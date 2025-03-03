import marimo

__generated_with = "0.11.13"
app = marimo.App(width="full")


@app.cell
def _():
    import pyfraglib as pf
    import pyfraglib.math as pm
    import marimo as mo
    import numpy as np
    import numpy.typing as npt

    from pyfraglib.math import fit_gmm, plot_gmm
    from scipy.optimize import LinearConstraint, minimize
    return LinearConstraint, fit_gmm, minimize, mo, np, npt, pf, plot_gmm, pm


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        The probability density function of two 1D Gaussians with variances $\sigma_1^2$, $\sigma_2^2$, means $\mu_1$, $\mu_2$, and mixture fraction $\pi$:

        \begin{equation}
            P(x) = \pi \cdot \mathcal{N}(x|\mu_1, \sigma_1^2) + (1-\pi) \cdot \mathcal{N}(x|\mu_2, \sigma_2^2)
        \end{equation}

        As a second piece, we need the negative log-likelihood:

        \begin{equation}
            \text{NLL} = -\sum{}{}{\log{P(x)}}
        \end{equation}

        Using pretty much any optimizer, we can then fit a GMM by minimizing the NLL:

        \begin{equation}
            \argmin_{\sigma_1^2, \sigma_2^2, \pi}{(\text{NLL})}
        \end{equation}
        """
    )
    return


@app.cell
def _(fit_gmm, np, npt, plot_gmm):
    _n: int = 10000
    mean1, mean2, mean3 = (160.0, 320.0, 520.0)
    std1, std2, std3 = (13.5, 20.0, 55.0)
    pi1, pi2, pi3 = (0.6, 0.3, 0.1)
    data1: npt.NDArray[np.float64] = np.random.normal(loc=mean1, scale=std1, size=int(_n * pi1))
    data2: npt.NDArray[np.float64] = np.random.normal(loc=mean2, scale=std2, size=int(_n * pi2))
    data3: npt.NDArray[np.float64] = np.random.normal(loc=mean3, scale=std3, size=int(_n * pi3))
    data: npt.NDArray[np.float64] = np.hstack([data1, data2, data3])

    _, _n, params, _ = fit_gmm(data, "../configs/gmm_3.json")
    plot_gmm(data, _n, params, ".", "test")
    return (
        data,
        data1,
        data2,
        data3,
        mean1,
        mean2,
        mean3,
        params,
        pi1,
        pi2,
        pi3,
        std1,
        std2,
        std3,
    )


if __name__ == "__main__":
    app.run()
