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

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from scipy.optimize import minimize
from scipy.stats import norm
from typing import Final

gmm_bounds_mu1: Final[tuple[float, float]] = (80, 200)
gmm_bounds_mu2: Final[tuple[float, float]] = (250, 420)


# @NOTE(ds): `params' is: mu_1, mu_2, std_1, std_2, pi.
def gaussian_mixture(
    params: tuple[float, float, float, float, float],
    data: list[float],
) -> npt.NDArray[np.float64]:
    m1, m2, std1, std2, pi = params

    pdf1 = norm.pdf(data, loc=m1, scale=std1)
    pdf2 = norm.pdf(data, loc=m2, scale=std2)

    pdf = pi * pdf1 + (1 - pi) * pdf2  # type: ignore
    return pdf  # type: ignore


# @NOTE(ds): `params' is: mu_1, mu_2, std_1, std_2, pi.
def negative_log_likelihood(
    params: tuple[float, float, float, float, float],
    data: list[float],
) -> npt.NDArray[np.float64]:
    pdf: npt.NDArray[np.float64] = gaussian_mixture(params, data)
    epsilon = 1e-10
    return -np.sum(np.log(pdf + epsilon))  # type: ignore


# @NOTE(ds): Fit a GMM of 2 1D Gaussians with bounds on the 5 free parameters.
# Right now, the bounds are hard-coded because we can make reasonable
# assumptions that should hold for (almost) all types of fragmentomics
# datasets.
def fit_gmm(data: npt.NDArray[np.float64]) -> list[float]:
    initial_params: list[float] = [
        167, 2*167,  # mu_1, mu_2
        5.0, 5.0,  # std_1, std_2
        0.5  # pi
    ]
    bounds: list[tuple[float, float]] = [
        gmm_bounds_mu1, gmm_bounds_mu2,
        (0.1, 50), (0.1, 100),  # std_1, std_2
        (0, 1)  # pi
    ]

    result = minimize(  # type: ignore
        negative_log_likelihood,
        initial_params,
        args=data,
        bounds=bounds,
        method='L-BFGS-B'
    )

    # Optimized parameters as [m1, m2, std1, std2, pi].
    return result.x  # type: ignore


# @NOTE(ds): Given 2 means and 2 standard deviations, we plot a histogram of
# `data' and overlay 2 normals. The parameters for the normals will most likely
# come from a GMM, that's why the function is named like this. `name' and
# `out_dir' are concatenated into a destination filepath.
def plot_gmm(data: npt.NDArray[np.float64], m1: float, m2: float, std1: float,
             std2: float, pi: float, out_dir: str, name: str) -> None:
    sample_size: int = len(data)
    sample_min: np.float64 = np.min(data)
    sample_max: np.float64 = np.max(data)
    bins: list[float] = list(
        np.arange(sample_min-30, sample_max+30)  # type: ignore
    )
    x: npt.NDArray[np.float64] = np.linspace(sample_min, sample_max, 1000)

    pdf1: npt.NDArray[np.float64] = \
        pi * norm.pdf(x, loc=m1, scale=std1)  # type: ignore
    pdf2: npt.NDArray[np.float64] = \
        (1 - pi) * norm.pdf(x, loc=m2, scale=std2)  # type: ignore
    pdf_gmm: npt.NDArray[np.float64] = pdf1 + pdf2

    fig = plt.figure()

    plt.hist(data, bins=bins, density=True, alpha=0.5, color="gray")
    plt.plot(x, pdf1, color="blue", linestyle="-.",
             label=r"$\mu_1={:.4}$, $\sigma_1={:.4}$".format(float(m1), std1))
    plt.plot(x, pdf2, color="orange", linestyle="--",
             label=r"$\mu_2={:.4}$, $\sigma_2={:.4}$".format(float(m2), std2))
    plt.plot(x, pdf_gmm, color="red", linewidth=2,
             label=r"GMM fit, $\pi={:.4}$".format(pi))
    plt.xlabel("Fragment Lengths")
    plt.ylabel("Density")
    plt.legend()
    plt.title("GMM for {}, n={}\n"
              r"$\mu_1 \in [{}, {}]$, $\mu_2 \in [{}, {}]$".format(
                  name, sample_size, gmm_bounds_mu1[0], gmm_bounds_mu1[1],
                  gmm_bounds_mu2[0], gmm_bounds_mu2[1]))

    outpath: str = \
        os.path.join(out_dir, "{}_gmm_frags_len.png".format(name))
    fig.savefig(outpath, dpi=fig.dpi)
