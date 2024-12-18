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

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from pyfraglib.core import fail
from scipy.optimize import LinearConstraint, minimize
from scipy.stats import norm


# @NOTE(ds): The ordering of `params' is as follows ('n' = # of Gaussians):
# mu_1, ..., mu_n, std_1, ... std_n, pi_1, ..., pi_n. The caller is responsible
# for ensuring that `sum(pi_i) == 1.0', and that `params' holds exactly n*3
# parameters.
def gaussian_mixture(
    params: list[float], n: int, data: list[float],
) -> npt.NDArray[np.float64]:
    assert n > 0, "Cannot mix <= 0 Gaussians."
    assert len(params) == 3*n, "Not enough mixture parameters."

    pdf: npt.NDArray[np.float64] | None = None
    for idx in range(n):
        m, s, p = params[idx], params[n+idx], params[2*n+idx]
        if pdf is None:
            pdf = norm.pdf(data, loc=m, scale=s) * p  # type: ignore
        else:
            pdf += norm.pdf(data, loc=m, scale=s) * p  # type: ignore
    return pdf  # type: ignore


# @NOTE(ds): `params' is as described above.
def negative_log_likelihood(
    params: list[float], n: int, data: list[float],
) -> npt.NDArray[np.float64]:
    # @NOTE(ds): We clip the mixture fractions to avoid numerical instability
    # issues.
    params[2*n:] = np.clip(params[2*n:], 1e-6, 1-1e-6)  # type: ignore
    pdf: npt.NDArray[np.float64] = gaussian_mixture(params, n, data)
    epsilon = 1e-10
    return -np.sum(np.log(pdf + epsilon))  # type: ignore


# @NOTE(ds): We do some basic validation, but don't exhaust possible error
# conditions in our sanity checks.
def read_gmm_config(
    config_filepath: str
) -> tuple[int, list[float], list[tuple[float, float]]]:
    with open(config_filepath, "r") as config_file:
        config: dict[str, object] = json.load(config_file)

        _n: object = config.get("number_of_gaussians")
        if not _n or type(_n) is not int or _n <= 0:
            fail("`number_of_gaussians' in GMM config must be an int > 0")
        n: int = int(_n)

        _ml: object = config.get("means_lower_bounds")
        _mu: object = config.get("means_upper_bounds")
        _sl: object = config.get("std_lower_bounds")
        _su: object = config.get("std_upper_bounds")
        _im: object = config.get("initial_means")
        _ip: object = config.get("initial_pis")

        if not _ml or not _mu or not _sl or not _su or not _im or not _ip:
            fail("missing required bounds array in GMM config")
        if (type(_ml) is not list or type(_mu) is not list or
                type(_sl) is not list or type(_su) is not list or
                type(_im) is not list or type(_ip) is not list):
            fail("bounds in GMM config must be float arrays")
        if (len(_ml) != n or len(_mu) != n or
                len(_sl) != n or len(_su) != n or
                len(_im) != n or len(_ip) != n):
            fail("bounds in GMM config must be float arrays")

        mean_bounds: list[tuple[float, float]] = []
        std_bounds: list[tuple[float, float]] = []
        mean_initial: list[float] = []
        std_initial: list[float] = []
        p_initial: list[float] = []
        var_zip: zip[tuple[float, float, float, float, float, float]] = \
            zip(_ml, _mu, _sl, _su, _im, _ip)
        for ml, mu, sl, su, im, ip in var_zip:
            mean_bounds.append((ml, mu))
            std_bounds.append((sl, su))
            mean_initial.append(im)
            p_initial.append(ip)

            # @NOTE(ds): Not configurable right now.
            std_initial.append((sl+su)/2)

        initials: list[float] = mean_initial + std_initial + p_initial
        p_bounds: list[tuple[float, float]] = [(0.0, 1.0)] * n
        bounds: list[tuple[float, float]] = mean_bounds + std_bounds + p_bounds

        return n, initials, bounds


def hessian(
    params: list[float], _n: int, _data: list[float],
) -> npt.NDArray[np.float64]:
    return np.zeros((len(params), len(params)))


# @NOTE(ds): Fit a GMM of `n' 1D Gaussians with bounds on the free parameters.
# The bounds are read from a configuration file. Optimized parameters are
# returned and ordered as follows:
# [m_1, m_2, ..., std_1, std_2, ..., pi_1, pi_2, ...].
def fit_gmm(
    data: npt.NDArray[np.float64], config_filepath: str
) -> tuple[object, int, list[float]]:
    initial_params: list[float]
    bounds: list[tuple[float, float]]
    n, initial_params, bounds = read_gmm_config(config_filepath)

    lin_constraint: LinearConstraint = LinearConstraint(
        A=[[0, 0, 0, 0, 0, 0, 1, 1, 1]],  lb=1, ub=1  # type: ignore
    )

    result = minimize(  # type: ignore
        negative_log_likelihood,
        initial_params,
        args=(n, data),
        bounds=bounds,
        constraints=(lin_constraint, ),
        method="trust-constr",
        hess=hessian,
        options={
            "maxiter": 1e6,
        }
    )

    return (result, n, result.x)  # type: ignore


# @NOTE(ds): Given `n' means and standard deviations, we plot a histogram of
# `data' and overlay `n' normals. The parameters for the normals will most
# likely come from a GMM, that's why the function is named like this. `name'
# and `out_dir' are concatenated into a destination filepath.
def plot_gmm(
    data: npt.NDArray[np.float64], num_gaussians: int, params: list[float],
    out_dir: str, name: str
) -> None:
    sample_size: int = len(data)
    sample_min: np.float64 = np.min(data)
    sample_max: np.float64 = np.max(data)
    bins: list[float] = list(
        np.arange(sample_min-30, sample_max+30)  # type: ignore
    )
    x: npt.NDArray[np.float64] = np.linspace(sample_min, sample_max, 1000)

    pdf_gmm: npt.NDArray[np.float64] | None = None
    n: int = num_gaussians

    fig = plt.figure()
    for idx in range(num_gaussians):
        m, s, p = params[idx], params[n+idx], params[2*n+idx]
        this_pdf: npt.NDArray[np.float64] = \
            norm.pdf(x, loc=m, scale=s) * p  # type: ignore
        plt.plot(x, this_pdf, linestyle="--",
                 label=r"$\mu_{}={:.4}$, $\sigma_{}={:.4}$, "
                       r"$\pi_{}={:.4}$".format(idx+1, m, idx+1, s, idx+1, p))
        if pdf_gmm is None:
            pdf_gmm = this_pdf
        else:
            pdf_gmm += this_pdf

    plt.hist(data, bins=bins, density=True, alpha=0.5, color="gray")
    plt.plot(x, pdf_gmm, color="red", linewidth=2, label="GMM fit")
    plt.xlabel("Fragment Lengths")
    plt.ylabel("Density")
    plt.legend()
    plt.title("GMM for {}, n={}".format(name, sample_size))

    outpath: str = \
        os.path.join(out_dir, "{}_gmm_frags_len.png".format(name))
    fig.savefig(outpath, dpi=fig.dpi)
