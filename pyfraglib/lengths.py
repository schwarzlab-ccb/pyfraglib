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
import logging

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from typing import Final
from pyfraglib import FragmentList, FragmentCollection


def fragment_length_histogram(
    fragments: FragmentList, out_dir: str, name: str
) -> None:
    logger: logging.Logger = logging.getLogger("pyfraglib")
    logger.info(
        "creating individual fragment length plots in `{}'".format(out_dir)
    )

    frag_lengths: list[int] = [
        frag.length for frag in fragments if not frag.is_bogus
    ]
    num_frags: Final = len(frag_lengths)

    fig: matplotlib.figure.Figure = plt.figure()
    outpath: str = os.path.join(out_dir, "{}_frags_len_histo.png".format(name))

    plt.hist(frag_lengths, bins=148)
    plt.title("{}: Fragment Length Distribution n={}".format(name, num_frags))
    plt.xlabel("Fragment Lengths")
    plt.ylabel("# of Fragments")

    fig.savefig(outpath, dpi=fig.dpi, transparent=True)

    frag_lengths_mut: list[int] = [
        frag.length for frag in fragments if
        (not frag.is_bogus and frag.is_mutated)
    ]
    frag_lengths_wt: list[int] = [
        frag.length for frag in fragments if
        (not frag.is_bogus and not frag.is_mutated)
    ]

    for do_density in (True, False):
        ptype: str = "density" if do_density else "histo"
        fig = plt.figure()
        outpath = os.path.join(
            out_dir, "{}_mut_frags_len_{}.png".format(name, ptype))

        plt.hist(frag_lengths_wt, bins=148, alpha=0.4, density=do_density,
                 label="wildtype fragments")
        plt.hist(frag_lengths_mut, bins=148, alpha=0.4, density=do_density,
                 label="mutated fragments")
        plt.title("{}: Fragment Length Distribution n={}".format(
            name, num_frags))
        plt.xlabel("Fragment Lengths")
        plt.ylabel("Density" if do_density else "# of Fragments")
        plt.legend()

        fig.savefig(outpath, dpi=fig.dpi, transparent=True)

    # Another option for plotting is a KDE plot.
    fig = plt.figure()
    outpath = os.path.join(out_dir, "{}_mut_frags_len_kde.png".format(name))

    sns.kdeplot(frag_lengths_wt, bw_adjust=1,
                label="wildtype, n={}".format(len(frag_lengths_wt)))
    sns.kdeplot(frag_lengths_mut, bw_adjust=1,
                label="mutated, n={}".format(len(frag_lengths_mut)))
    plt.title("{}: Fragment Length Distribution n={}".format(name, num_frags))
    plt.xlabel("Fragment Lengths")
    plt.ylabel("Density")
    plt.legend()

    fig.savefig(outpath, dpi=fig.dpi, transparent=True)


def fragment_length_histograms(
    fragments: FragmentCollection, out_dir: str
) -> None:
    logger: logging.Logger = logging.getLogger("pyfraglib")
    logger.info(
        "creating fragment length plots for all samples in `{}'".format(
            out_dir))

    name: str
    frag_list: FragmentList
    outpath: str = os.path.join(out_dir, "all_bams_frags_len_histo.png")
    fig: matplotlib.figure.Figure = plt.figure()

    for name, frag_list in fragments:
        frag_lengths: list[int] = [
            frag.length for frag in frag_list if not frag.is_bogus
        ]
        num_frags: int = len(frag_lengths)
        sns.kdeplot(frag_lengths, bw_adjust=1, alpha=0.8,
                    label="{} (n={})".format(name, num_frags))

    plt.title("Fragment Length Distribution")
    plt.xlabel("Fragment Lengths")
    plt.ylabel("# of Fragments")
    plt.legend()

    fig.savefig(outpath, dpi=fig.dpi, transparent=True)
