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
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize
from scipy.stats import norm


def gaussian_mixture(
    params: tuple[float, float, float],
    data: list[float],
    fixed_mus: tuple[float, float]
) -> npt.NDArray[np.float64]:
    pi, std1, std2 = params
    m1, m2 = fixed_mus

    pdf1 = norm.pdf(data, loc=m1, scale=std1)
    pdf2 = norm.pdf(data, loc=m2, scale=std2)

    pdf = pi * pdf1 + (1 - pi) * pdf2  # type: ignore
    return pdf  # type: ignore


def negative_log_likelihood(
    params: tuple[float, float, float],
    data: list[float],
    fixed_mus: tuple[float, float]
) -> npt.NDArray[np.float64]:
    pdf: npt.NDArray[np.float64] = gaussian_mixture(params, data, fixed_mus)
    epsilon = 1e-10
    return -np.sum(np.log(pdf + epsilon))  # type: ignore


def fit_gmm(
    m1: float, m2: float,
    data: npt.NDArray[np.float64]
) -> list[float]:
    fixed_means: list[float] = [m1, m2]

    # `initial_params' and `bounds' must be in order "pi", "std1", "std2".
    initial_params: list[float] = [0.5, 1.0, 1.0]
    bounds: list[tuple[float, float]] = [(0, 1), (1e-3, 1000), (1e-3, 1000)]

    result = minimize(  # type: ignore
        negative_log_likelihood,
        initial_params,
        args=(data, fixed_means),
        bounds=bounds,
        method='L-BFGS-B'
    )

    # Optimized parameters as [pi, std1, std2].
    return result.x  # type: ignore


# @NOTE(ds): Given 2 means and 2 standard deviations, we plot a histogram of
# `data' and overlay 2 normals. The parameters for the normals will most likely
# come from a GMM, that's why the function is named like this.
def plot_gmm(data: npt.NDArray[np.float64], m1: float, m2: float,
             pi: float, std1: float, std2: float, bins: int = 50) -> None:
    x: npt.NDArray[np.float64] = np.linspace(np.min(data), np.max(data), 1000)

    pdf1: npt.NDArray[np.float64] = \
        pi * norm.pdf(x, loc=m1, scale=std1)  # type: ignore
    pdf2: npt.NDArray[np.float64] = \
        (1 - pi) * norm.pdf(x, loc=m2, scale=std2)  # type: ignore
    pdf_gmm: npt.NDArray[np.float64] = pdf1 + pdf2

    plt.hist(data, bins=bins, density=True, alpha=0.5, color="gray")
    plt.plot(x, pdf1, label="Gaussian 1", color="blue", linestyle="-.")
    plt.plot(x, pdf2, label="Gaussian 2", color="orange", linestyle="--")
    plt.plot(x, pdf_gmm, label="GMM (fitted)", color="red", linewidth=2)
    plt.xlabel("Data value")
    plt.ylabel("Density")
    plt.legend()
    plt.title("Gaussian Mixture Model Fit")
    plt.show()  # type: ignore
