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
import os
import json
import logging
import itertools

import matplotlib
import matplotlib.pyplot as plt

from collections import defaultdict
from pyfraglib.core import fail
from pyfraglib.fragment import FragmentList
from pyfraglib import get_logger


def fragments_per_chromosome_barplot(
    fragments: FragmentList, out_dir: str, name: str
) -> None:
    logger: logging.Logger = get_logger()
    logger.info(
        "creating fragment-per-chromosome histograms in `{}'".format(out_dir)
    )

    chrom_map_mut: defaultdict[str, int] = defaultdict(int)
    chrom_map_wt:  defaultdict[str, int] = defaultdict(int)
    num_frags = 0
    for frag in fragments:
        if frag.is_bogus:
            continue

        num_frags += 1
        if frag.is_mutated:
            chrom_map_mut[frag.chrom] += 1
            _ = chrom_map_wt[frag.chrom]
        else:
            chrom_map_wt[frag.chrom] += 1
            _ = chrom_map_mut[frag.chrom]

    fig = plt.figure()
    outpath = os.path.join(out_dir, "{}_mut_frags_per_chrom.png".format(name))

    assert chrom_map_wt.keys() == chrom_map_mut.keys()
    plt.bar(list(chrom_map_wt.keys()), list(chrom_map_wt.values()),
            label="wildtype fragments")
    if sum(chrom_map_mut.values()) != 0:
        plt.bar(list(chrom_map_mut.keys()), list(chrom_map_mut.values()),
                bottom=list(chrom_map_wt.values()), label="mutated fragments")

    plt.xticks(rotation=45)  # better rotate bc. we do not homogenize chr names
    plt.title("{}: Fragments per chromosome n={}".format(name, num_frags))
    plt.xlabel("Chromosomes")
    plt.ylabel("# of Fragments")
    plt.legend()

    fig.savefig(outpath, dpi=900)


def end_motifs_barplot(
    fragments: FragmentList, out_dir: str, name: str, kmer_len: int
) -> None:
    logger: logging.Logger = get_logger()
    logger.info(
        "creating end motif barplots in `{}'".format(out_dir)
    )

    motifs_5p: defaultdict[str, int]
    motifs_3p: defaultdict[str, int]
    num_frags: int

    motifs_5p, motifs_3p, num_frags = fragments.count_endmotifs(kmer_len)
    assert motifs_5p.keys() == motifs_3p.keys()

    for label, motifs in [("5", motifs_5p), ("3", motifs_3p)]:
        mzip: zip[tuple[str, int]] = zip(motifs.keys(), motifs.values())
        szip: list[tuple[str, int]] = sorted(mzip, reverse=False)

        motifs_names: list[str] = []
        motifs_vals: list[int] = []
        m_name: str
        m_val: int
        for m_name, m_val in szip:
            motifs_names.append(m_name)
            motifs_vals.append(m_val)

        # @NOTE(ds): We assume sections of equal lengths!
        max_count: int = max(motifs_vals)
        sec_len: int = len(motifs_names) // 4

        fig: matplotlib.figure.Figure = plt.figure()
        outpath: str = os.path.join(
            out_dir, "{}_k{}_{}p_frag_end_motifs.png".format(
                name, kmer_len, label)
        )

        colors: list[str] = ["red", "green", "blue", "gold"]
        bases: list[str] = ["A", "C", "G", "T"]
        colors_bars: list[str] = list(itertools.chain(
            *[sec_len*[color] for color in colors]))

        plt.bar(motifs_names, motifs_vals, color=colors_bars)
        for num, (color, base) in enumerate(zip(colors, bases)):
            x_start: float = sec_len*num-0.5
            x_end: float = sec_len*(num+1)-0.5
            x_mid: float = (x_start + x_end) / 2
            plt.hlines(y=max_count*1.15, color=color, linewidth=3,
                       xmin=x_start, xmax=x_end)
            plt.text(x_mid, max_count*1.12, base, ha="center", va="top",
                     color=color, fontsize=18)

        plt.title("{}: {}' End motif distribution n={}".format(
            name, label, num_frags))
        plt.xlabel("{}-mer End Motifs".format(kmer_len))
        plt.ylabel("# of Motifs")
        plt.xticks(rotation=80, fontsize=0.6)
        plt.tick_params(axis="x", length=0, pad=0.4)

        fig.savefig(outpath, dpi=1200, bbox_inches="tight")


def log_stats(
    fragments: FragmentList, logger: logging.Logger,
    out_dir: str, name: str
) -> None:
    num_bogus_frags: int = fragments.count_bogus_fragments()
    num_mut_frags: int = fragments.count_mutated_fragments()
    num_frags: int = fragments.length()

    if num_frags == 0:
        fail("no fragments found")

    logger.info(
        "identified {} fragments (bogus: {} ({}%), mutated: {} ({}%))".format(
            num_frags, num_bogus_frags,
            round(100*num_bogus_frags/num_frags, 3),
            num_mut_frags,
            round(100*num_mut_frags/num_frags, 3))
    )

    outpath: str = os.path.join(out_dir, "{}_frag_stats.json".format(name))
    with open(outpath, "w") as file:
        data: dict[str, object] = {
            "number_of_fragments": num_frags,
            "number_of_bogus_fragments": num_bogus_frags,
            "number_of_mutated_fragments": num_mut_frags
        }
        json.dump(data, file)
