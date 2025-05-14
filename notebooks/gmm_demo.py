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
            P(x) = \pi_1 \cdot \mathcal{N}(x|\mu_1, \sigma_1) + \pi_2 \cdot \mathcal{N}(x|\mu_2, \sigma_2) + \pi_3 \cdot \mathcal{N}(x|\mu_3, \sigma_3)
        \end{equation}

        As a second piece, we need the negative log-likelihood:

        \begin{equation}
            \text{NLL}(\mathcal{X}) = -\sum_{x \in \mathcal{X}}^{}{\log{P(x)}}
        \end{equation}

        Using pretty much any optimizer, we can then fit a GMM by minimizing the NLL:

        \begin{equation}
            \argmin_{\mu_1, \mu_2, \mu_3, \sigma_1, \sigma_2, \sigma_3, \pi_1, \pi_2, \pi_3}{(\text{NLL})},
        \end{equation}

        Such that:

        \begin{equation}
            \sum_{i=1, \pi_i \ge 0}^{3}{pi_i} = 1
        \end{equation}

        \begin{equation}
            \text{for } i \in \{1, 2, 3\}: \sigma_i \in (0, 50]; \mu_1 \in [50, 200]; \mu_2 \in [200, 350]; \mu_3 \in [350, 600]
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
    plot_gmm(data, _n, params, "/home/daniel/", "test")
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
