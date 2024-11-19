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
                           get_chromosome_length, hg19_chromosomes
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
    fragments: FragmentList, regions: pysam.TabixFile, win_size: int = 120
) -> pd.DataFrame:
    return windowed_protection_score_fast(fragments, regions, win_size)


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
    fragments: FragmentList, regions: pysam.TabixFile, win_size: int = 120
) -> pd.DataFrame:
    assert win_size > 0

    chromosome_map: dict[str, npt.NDArray[np.int64]] = create_chromosome_map()

    for fragment in fragments:
        if fragment.is_bogus:
            continue

        frag_start: int = fragment.start_pos
        frag_end: int = fragment.end_pos  # 1 past end
        frag_len: int = fragment.length
        frag_chrom: str = fragment.chrom

        win_half: int = win_size // 2
        istart: int = max(frag_start - win_half, 0)
        iend: int = min(frag_end + win_half, get_chromosome_length(frag_chrom))

        if frag_len < win_size:
            chromosome_map[frag_chrom][istart:iend] -= 1  # type: ignore
        else:
            ostart: int = frag_start + win_half
            oend: int = frag_end - win_half  # 1 past end

            chromosome_map[frag_chrom][istart:ostart] -= 1  # type: ignore
            chromosome_map[frag_chrom][ostart:oend] += 1  # type: ignore
            chromosome_map[frag_chrom][oend:iend] -= 1  # type: ignore

    return chromosome_map_to_df(chromosome_map, regions)


# @NOTE(ds): We assume a correctly formatted BED file to be provided for the
# `region' argument. An `info' field can be used to carry additional data over
# the the resulting dataframe (e.g. gene names).
def windowed_protection_score_slow(
    fragments: FragmentList, regions: pysam.TabixFile,
    win_size: int = 120, step_size: int = 1,
) -> pd.DataFrame:
    assert win_size > 0
    assert step_size > 0

    col_names: list[str] = ["chrom", "pos", "abs_pos", "wps", "info"]
    wps_df: pd.DataFrame = pd.DataFrame(
        None, index=range(precalc_size(regions, step_size)), columns=col_names
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
        for win_start in range(int(istart)-60, int(iend)-60+1, step_size):
            win_end: int = win_start + win_size  # right-exclusive
            win_mid: int = win_start + win_size // 2

            if not cur_chrom:
                cur_chrom = chrom
            elif chrom != cur_chrom:
                cum_pos += get_chromosome_length(cur_chrom)
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


def precalc_size(regions: pysam.TabixFile, step_size: int) -> int:
    it: int = 0
    for region in regions.fetch():  # type: ignore
        istart: str
        iend: str
        _, istart, iend, _ = region.split()  # type: ignore
        it += len(range(int(istart), int(iend)+1, step_size))
    regions.reset()
    return it


def create_chromosome_map() -> dict[str, npt.NDArray[np.int64]]:
    chrom_map: dict[str, npt.NDArray[np.int64]] = dict()
    name: str
    length: int

    for name, length, _, _ in hg19_chromosomes:
        chrom_map[name] = np.zeros(length, dtype=np.int64)

    return chrom_map


# @NOTE(ds): The input BED file must be sorted with respect to chromosomes.
# E.g. listing chr1 regions, chr3 regions, and chr2 regions in this order will
# produce incorrect results with respect to the absolute genome position!
def chromosome_map_to_df(
    chrom_map: dict[str, npt.NDArray[np.int64]], regions: pysam.TabixFile
) -> pd.DataFrame:
    col_names: list[str] = ["chrom", "pos", "abs_pos", "wps", "info"]
    wps_df: pd.DataFrame = pd.DataFrame(
        None, index=range(precalc_size(regions, 1)),
        columns=col_names
    )

    region: str
    it: int = 0
    cur_chrom: str | None = None
    cum_pos: int = 0
    chromosome: npt.NDArray[np.int64] | None = None

    for region in regions.fetch():  # type: ignore
        chrom: str
        istart: str
        iend: str
        info: str

        chrom, istart, iend, info = region.split()

        if not cur_chrom or chromosome is None:
            cur_chrom = chrom
            chromosome = chrom_map[chrom]
        elif chrom != cur_chrom:
            cum_pos += get_chromosome_length(cur_chrom)
            cur_chrom = chrom
            chromosome = chrom_map[chrom]

        rlen: int = int(iend) - int(istart)  # region length - 1
        rel_pos: npt.NDArray[np.int64] = np.arange(int(istart), int(iend) + 1)

        wps_df.loc[it:(it+rlen), "chrom"] = chrom
        wps_df.loc[it:(it+rlen), "pos"] = rel_pos
        wps_df.loc[it:(it+rlen), "abs_pos"] = rel_pos + cum_pos  # type: ignore
        wps_df.loc[it:(it+rlen), "wps"] = chromosome[rel_pos]  # type: ignore
        wps_df.loc[it:(it+rlen), "info"] = info

        it += rlen + 1

    regions.reset()

    return wps_df


def wps_scatter_plot(df: pd.DataFrame, name: str, out_dir: str) -> None:
    outpath: str = os.path.join(out_dir, "{}_wps.png".format(name))
    logging.getLogger("pyfraglib").info(
        "saving WPS plot for {} to `{}'".format(name, outpath)
    )

    fig: matplotlib.figure.Figure = plt.figure()

    plt.plot(df["abs_pos"], df["wps"])  # type: ignore
    plt.title("Windowed Protection Score for {}".format(name))
    plt.xlabel("Genomic Position")
    plt.ylabel("WPS")

    fig.savefig(outpath, dpi=900)
