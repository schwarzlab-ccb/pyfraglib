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

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import seaborn as sns

from typing import Final
from pyfraglib import FragmentList
from pyfraglib.math import fit_gmm, plot_gmm


def fragment_length_plot(
    fragments: FragmentList, out_dir: str, name: str
) -> None:
    logger: logging.Logger = logging.getLogger("pyfraglib")
    logger.info(
        "creating fragment length plots in `{}'".format(out_dir)
    )

    frag_lengths_mut: list[int] = []
    frag_lengths_wt: list[int] = []

    for frag in fragments:
        if frag.is_bogus:
            continue
        if frag.is_mutated:
            frag_lengths_mut.append(frag.length)
        else:
            frag_lengths_wt.append(frag.length)

    num_frags: Final[int] = len(frag_lengths_mut) + len(frag_lengths_wt)

    fig = plt.figure()
    outpath: str = \
        os.path.join(out_dir, "{}_mut_frags_len_kde.png".format(name))

    sns.kdeplot(frag_lengths_wt, bw_adjust=1,
                label="wildtype, n={}".format(len(frag_lengths_wt)))
    sns.kdeplot(frag_lengths_mut, bw_adjust=1,
                label="mutated, n={}".format(len(frag_lengths_mut)))
    plt.title("{}: Fragment Length Distribution n={}".format(name, num_frags))
    plt.xlabel("Fragment Lengths")
    plt.ylabel("Density")
    plt.legend()

    fig.savefig(outpath, dpi=fig.dpi)


def fragment_length_gmm(fragments: FragmentList,
                        out_dir: str, name: str) -> None:
    logger: logging.Logger = logging.getLogger("pyfraglib")
    logger.info(
        "fitting GMM, writing results to `{}'".format(out_dir)
    )

    frag_lens: npt.NDArray[np.float64] = np.array(
        [frag.length for frag in fragments if not frag.is_bogus]
    )

    m1, m2 = 167, 2*167
    [pi, std1, std2] = fit_gmm(m1, m2, frag_lens)
    plot_gmm(frag_lens, m1, m2, pi, std1, std2, out_dir, name)
