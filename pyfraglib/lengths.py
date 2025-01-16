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

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import seaborn as sns

from typing import Final
from pyfraglib import FragmentList, fail
from pyfraglib.math import fit_gmm, plot_gmm, goodness_of_fit_stats


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
    plt.title("{}: Length Distribution n={}".format(name, num_frags))
    plt.xlabel("Fragment Lengths")
    plt.ylabel("Density")
    plt.legend()

    fig.savefig(outpath, dpi=fig.dpi)


def fragment_length_gmm(fragments: FragmentList, config_filepath: str,
                        out_dir: str, name: str) -> None:
    logger: logging.Logger = logging.getLogger("pyfraglib")
    logger.info(
        "fitting GMM based on config file `{}', writing results to `{}'"
        .format(config_filepath, out_dir)
    )

    frag_lens: npt.NDArray[np.float64] = np.array(
        [frag.length for frag in fragments if not frag.is_bogus]
    )

    try:
        opt_result, n, params, data = fit_gmm(frag_lens, config_filepath)
    except Exception:
        fail("fitting the GMM failed, probably due to non-convergence")

    if opt_result.success:  # type: ignore
        logger.info("successfully fitted GMM")
    else:
        fail("fitting the GMM failed due to non-convergence")

    plot_gmm(data, n, params, out_dir, name)

    gof: dict[str, object] = goodness_of_fit_stats(data, params, n)
    write_gmm_params(n, params,
                     opt_result.fun,  opt_result.success,  # type: ignore
                     gof, out_dir, name)


def write_gmm_params(
    n: int, params: list[float], obj_val: float, conv: bool,
    goodness_of_fit: dict[str, object], out_dir: str, name: str
) -> None:
    outpath: str = \
        os.path.join(out_dir, "{}_gmm_frags_len.json".format(name))
    with open(outpath, "w") as file:
        data: dict[str, object] = {
            "number_of_gaussians": n,
            "objective_value": obj_val,
            "converged": conv,
            "estimated_means": list(params[:n]),
            "estimated_stds": list(params[n:2*n]),
            "estimated_pis": list(params[2*n:])
        }
        data |= goodness_of_fit
        json.dump(data, file)
