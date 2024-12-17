# This file is part of `pyfraglib`, a software suite to calculate fragmentomics
# features from cfDNA and perform downstream analyses.
#
# Copyright (C) 2024 Daniel Schütte, daniel.schuette@iccb-cologne.org
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details. You should have received a copy of the GNU General Public
# License along with this program. If not, see <https://www.gnu.org/licenses/>.
import logging
import typing
import os
import pysam
import matplotlib

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from collections import defaultdict
from intervaltree import Interval  # type: ignore
from pyfraglib.core import shannon_entropy, simpson_index, fail, \
                           get_chromosome_length, hg19_chromosomes, \
                           hg38_chromosomes, homogenize_contig_name
from pyfraglib.fragment import FragmentList, IntervalTable


# NOTE(ds): Return the MDS for 5' and 3' fragments.
def motif_diversity(
    fragments: FragmentList, name: str, index: str = "shannon"
) -> tuple[float, float]:
    logger: logging.Logger = logging.getLogger("pyfraglib")

    index_func: typing.Callable[[list[float]], float]
    if index == "shannon":
        index_func = shannon_entropy
    elif index == "simpson":
        index_func = simpson_index
    else:
        fail("motif_diversity(): unknown index function `{}'".format(index))

    logger.info(
        "calculating motif diversity score ({}) for {} ".format(index, name)
    )

    ends_5p: defaultdict[str, int] = defaultdict(int)
    ends_3p: defaultdict[str, int] = defaultdict(int)

    ends_5p, ends_3p, num_frags = fragments.count_endmotifs(kmer_len=4)

    total_5p: int = sum(ends_5p.values())
    proportions_5p: list[float] = [val / total_5p for val in ends_5p.values()]
    index_5p: float = index_func(proportions_5p)

    total_3p: int = sum(ends_3p.values())
    proportions_3p: list[float] = [val / total_3p for val in ends_3p.values()]
    index_3p: float = index_func(proportions_3p)

    return (index_5p, index_3p)


# @NOTE(ds): This is a simply dispatch function to facilitate selection of a
# WPS implementation.
def windowed_protection_score(
    fragments: FragmentList, regions: pysam.TabixFile, win_size: int = 120,
    genome: str = "hg19"
) -> pd.DataFrame:
    return windowed_protection_score_fast(fragments, regions, win_size, genome)


# @NOTE(ds): We removed the `step_size' argument because this algorithm does
# not need that. It iterates over fragments instead of genomic regions:
#
#  > case 1 (fragment_size < window_size):
#       [frag_start - win_size/2, frag_end + win_size/2] <- -1
#
# > case 2 (fragment_size >= window_size):
#       [frag_start - win_size/2, frag_start + win_size/2) <- -1
#       [frag_start + win_size/2, frag_end   - win_size/2] <- +1
#       (frag_end   - win_size/2, frag_end   + win_size/2] <- -1
#
# The implementation below tries to be extremely careful about interval
# definitions to not introduce 1-off errors.
def windowed_protection_score_fast(
    fragments: FragmentList, regions: pysam.TabixFile, win_size: int = 120,
    genome: str = "hg19"
) -> pd.DataFrame:
    assert win_size > 0

    chromosome_map_wps: dict[str, npt.NDArray[np.int64]] = \
        create_chromosome_map(genome)
    chromosome_map_depth: dict[str, npt.NDArray[np.int64]] = \
        create_chromosome_map(genome)

    for fragment in fragments:
        if fragment.is_bogus:
            continue

        frag_start: int = fragment.start_pos
        frag_end: int = fragment.end_pos  # 1 past end
        frag_len: int = fragment.length
        frag_chrom: str = homogenize_contig_name(fragment.chrom)

        win_half: int = win_size // 2
        istart: int = max(frag_start - win_half, 0)
        iend: int = min(frag_end + win_half,
                        get_chromosome_length(frag_chrom, genome))

        chromosome_map_depth[frag_chrom][istart:iend] += 1  # type: ignore
        if frag_len < win_size:
            chromosome_map_wps[frag_chrom][istart:iend] -= 1  # type: ignore
        else:
            ostart: int = frag_start + win_half
            oend: int = frag_end - win_half  # 1 past end

            chromosome_map_wps[frag_chrom][istart:ostart] -= 1  # type: ignore
            chromosome_map_wps[frag_chrom][ostart:oend] += 1  # type: ignore
            chromosome_map_wps[frag_chrom][oend:iend] -= 1  # type: ignore

    return chromosome_maps_to_df(chromosome_map_wps, chromosome_map_depth,
                                 regions, genome)


# @DEPRECATED(ds): This function should no longer be used!
#
# @NOTE(ds): We assume a correctly formatted BED file to be provided for the
# `region' argument. An `info' field can be used to carry additional data over
# the the resulting dataframe (e.g. gene names).
def windowed_protection_score_slow(
    fragments: FragmentList, regions: pysam.TabixFile,
    win_size: int = 120, genome: str = "hg19"
) -> pd.DataFrame:
    assert win_size > 0

    col_names: list[str] = ["chrom", "pos", "abs_pos", "wps", "info"]
    wps_df: pd.DataFrame = pd.DataFrame(
        None, index=range(precalc_size(regions)), columns=col_names
    )

    region: str
    interval_table: IntervalTable = fragments.to_interval_table()
    it: int = 0
    cum_pos: int = 0
    cur_chrom: str | None = None

    for region in regions.fetch():  # type: ignore
        chrom: str
        istart: str
        iend: str
        info: str

        chrom, istart, iend, info = region.split()
        logging.getLogger("pyfraglib").debug("WPS at {}".format(info))

        # @NOTE(ds): To calculate our scores over the same windows as we do
        # with the fast algorithm, we need to slightly adjust the windows.
        for win_start in range(int(istart)-60, int(iend)-60+1):
            win_end: int = win_start + win_size  # right-exclusive
            win_mid: int = win_start + win_size // 2

            if not cur_chrom:
                cur_chrom = chrom
            elif chrom != cur_chrom:
                cum_pos += get_chromosome_length(cur_chrom, genome)
                cur_chrom = chrom

            intervals: list[Interval]  # type: ignore
            intervals = interval_table.get_overlaps(
                chrom, win_start, win_end
            )

            num_spanning_frags: int = 0
            num_intersecting_frags: int = 0
            interval: Interval  # type: ignore
            for interval in intervals:  # type: ignore
                s: int = interval.begin  # type: ignore
                e: int = interval.end  # type: ignore
                if s > win_start or e < win_end:
                    num_intersecting_frags += 1
                else:
                    num_spanning_frags += 1

            new_row: list[str | float] = [
                chrom, win_mid, win_mid + cum_pos,
                num_spanning_frags - num_intersecting_frags,
                info
            ]
            wps_df.loc[it] = new_row
            it += 1

    regions.reset()

    logging.getLogger("pyfraglib").debug(
        "calculated WPS for {} positions".format(it)
    )
    return wps_df


def precalc_size(regions: pysam.TabixFile) -> int:
    it: int = 0
    for region in regions.fetch():  # type: ignore
        istart: str
        iend: str
        _, istart, iend, _ = region.split()  # type: ignore
        it += int(iend) - int(istart) + 1
    regions.reset()
    return it


def create_chromosome_map(
    genome: str = "hg19"
) -> dict[str, npt.NDArray[np.int64]]:
    chrom_map: dict[str, npt.NDArray[np.int64]] = dict()
    name: str
    length: int
    chromosomes: list[tuple[str, int, str, str]]

    if genome == "hg19":
        chromosomes = hg19_chromosomes
    elif genome == "hg38":
        chromosomes = hg38_chromosomes
    else:
        fail("unknown genome `{}' requested".format(genome))

    for name, length, _, _ in chromosomes:
        chrom_map[name] = np.zeros(length, dtype=np.int64)

    return chrom_map


# @NOTE(ds): The input BED file must be sorted with respect to chromosomes.
# E.g. listing chr1 regions, chr3 regions, and chr2 regions in this order will
# produce incorrect results with respect to the absolute genome position!
def chromosome_maps_to_df(
    chrom_map_wps: dict[str, npt.NDArray[np.int64]],
    chrom_map_depth: dict[str, npt.NDArray[np.int64]],
    regions: pysam.TabixFile, genome: str = "hg19"
) -> pd.DataFrame:
    col_names: list[str] = ["chrom", "pos", "abs_pos", "wps", "depth", "info"]
    output_df: pd.DataFrame = pd.DataFrame(
        None, index=range(precalc_size(regions)),
        columns=col_names
    )

    region: str
    it: int = 0
    cur_chrom: str | None = None
    cum_pos: int = 0
    chromosome_wps: npt.NDArray[np.int64] | None = None
    chromosome_depth: npt.NDArray[np.int64] | None = None

    for region in regions.fetch():  # type: ignore
        chrom: str
        istart: str
        iend: str
        info: str

        chrom, istart, iend, info = region.split()
        chrom = homogenize_contig_name(chrom)

        if not cur_chrom or chromosome_wps is None or chromosome_depth is None:
            cur_chrom = chrom
            chromosome_wps = chrom_map_wps[chrom]
            chromosome_depth = chrom_map_depth[chrom]
        elif chrom != cur_chrom:
            cum_pos += get_chromosome_length(cur_chrom, genome)
            cur_chrom = chrom
            chromosome_wps = chrom_map_wps[chrom]
            chromosome_depth = chrom_map_depth[chrom]

        rlen: int = int(iend) - int(istart)  # region length - 1
        rel_pos: npt.NDArray[np.int64] = np.arange(int(istart), int(iend) + 1)

        output_df.loc[it:(it+rlen), "chrom"] = chrom
        output_df.loc[it:(it+rlen), "pos"] = rel_pos
        output_df.loc[it:(it+rlen), "abs_pos"] = \
            rel_pos + cum_pos  # type: ignore
        output_df.loc[it:(it+rlen), "wps"] = \
            chromosome_wps[rel_pos]  # type: ignore
        output_df.loc[it:(it+rlen), "depth"] = \
            chromosome_depth[rel_pos]  # type: ignore
        output_df.loc[it:(it+rlen), "info"] = info

        it += rlen + 1

    regions.reset()

    return output_df


# @NOTE(ds): `exclude_chroms' must be in "1" format (not "chr1"). `region` is
# only applied if a single chromosome is selected (by excluding all chromosomes
# but one, that is).
def score_line_plot(
    df: pd.DataFrame, name: str, out_dir: str, score: str = "wps",
    exclude_chroms: list[str] = ["Y", "M"], region: tuple[int, int] = (0, 0),
    genome: str = "hg19"
) -> None:
    if score not in ["wps", "depth"]:
        fail("unknown score `{}' requested".format(score))

    outpath: str = os.path.join(out_dir, "{}_{}.png".format(name, score))
    logging.getLogger("pyfraglib").info(
        "saving score plot for {} to `{}'".format(name, outpath)
    )

    chromosomes: list[tuple[str, int, str, str]]
    if genome == "hg19":
        chromosomes = hg19_chromosomes
    elif genome == "hg38":
        chromosomes = hg38_chromosomes
    else:
        fail("unknown genome `{}' requested".format(genome))

    fig: matplotlib.figure.Figure = plt.figure()
    plt.title("Sample {}".format(name))
    plt.xlabel("Chromosomes")
    plt.ylabel(score.capitalize())

    cumsum: int = 0
    chrom_midpoints: list[int] = []
    chrom_names: list[str] = []
    color: str = "lightgrey"
    for name, length, _, _ in chromosomes:
        if name in exclude_chroms:
            continue

        chrom_midpoints.append(cumsum + length//2)
        chrom_names.append(name)

        plt.axvspan(cumsum, cumsum+length, color=color, alpha=0.3)
        if color == "lightgrey":
            color = "white"
        else:
            color = "lightgrey"

        condition = df["chrom"] == name  # type: ignore
        plt.plot(
            df.loc[condition, "pos"] + cumsum,  # type: ignore
            df.loc[condition, score],  # type: ignore
            color="#1f77b4"
        )

        cumsum += length

    if len(chrom_names) == 1:
        plt.xlim(region)
        plt.xticks([region[0], region[1], np.mean(region)],  # type: ignore
                   [region[0], region[1], chrom_names[0]])  # type: ignore
    else:
        plt.xlim(-cumsum*0.05, cumsum*1.05)
        staggered_chrom_names: list[str] = [
            f"{label}\n" if i % 2 == 0 else f"\n{label}"
            for i, label in enumerate(chrom_names)
        ]
        plt.xticks(chrom_midpoints, staggered_chrom_names, rotation=0)

    plt.tight_layout()
    fig.savefig(outpath, dpi=900)
