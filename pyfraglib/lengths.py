"""
Fragment length analysis and visualization
==========================================

This module provides tools for analyzing and visualizing fragment length
distributions. It includes functions for creating kernel density plots to
compare wildtype and mutated fragments, and for fitting Gaussian mixture models
(GMMs) to fragment length data.

The module supports:

- **Length distribution visualization**: Kernel density estimation plots
                                         comparing wildtype vs mutated fragment
                                         length distributions
- **Gaussian mixture modeling**: Automated fitting of multi-component normal
                                 distributions to fragment length data with
                                 convergence checking
- **Goodness-of-fit statistics**: Statistical evaluation of GMM fits including
                                  Kolmogorov-Smirnov tests and Jensen-Shannon
                                  divergence
- **Parameter export**: JSON serialization of fitted GMM parameters and
                        statistics

Key Features:
-------------

* Automated GMM fitting with configurable number of components
* Statistical validation of model fits
* Visualization of results

Example:
    Basic fragment length analysis workflow::

        from pyfraglib import Fragment, fragment_length_plot, \
                              fragment_length_gmm

        # Load fragments from BAM file
        fragments = Fragment.from_bam("sample.bam", "variants.vcf")

        # Create length distribution plots
        fragment_length_plot(fragments, "output/", "sample")

        # Fit Gaussian mixture model
        fragment_length_gmm(
            fragments, "configs/gmm_3.json", "output/", "sample"
        )

The analysis generates multiple outputs including visualizations, fitted
parameters, and goodness-of-fit statistics for downstream analysis and quality
control.

See Also
---------
    * :mod:`pyfraglib.math` - Core mathematical functions for GMM fitting
    * :mod:`pyfraglib.stats` - General statistical analysis functions
    * :doc:`../examples/index` - Complete analysis examples

License
-------
This file is part of ``pyfraglib``, a software suite to calculate
fragmentomics features from cfDNA and perform downstream analyses.

Copyright (C) 2025 Daniel Schütte, daniel.schuette@iccb-cologne.org

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details. You should have received a copy of the GNU General Public
License along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
import os
import json
import logging

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import seaborn as sns

from typing import Final
from pyfraglib import FragmentList, fail, get_logger
from pyfraglib.math import fit_gmm, plot_gmm, goodness_of_fit_stats


def fragment_length_plot(
    fragments: FragmentList, out_dir: str, name: str
) -> None:
    """
    Create kernel density estimation plots for fragment length distributions.

    Generates exploratory plots comparing fragment length distributions
    between wildtype and mutated fragments using kernel density estimation.
    The plot helps visualize differences in fragmentation patterns between
    normal and variant-carrying cfDNA fragments.

    Args:
        fragments: List of fragments to analyze. Should contain both
            wildtype and mutated fragments for meaningful comparison.
        out_dir: Output directory path where the plot will be saved.
            Directory will be created if it doesn't exist.
        name: Base name for the output file. The plot will be saved as
            ``{name}_frags_len_kde.png``.

    Returns:
        None. The function saves a PNG plot to the specified output directory.

    Note:
        * Bogus fragments (quality-filtered) are automatically excluded
        * Includes sample sizes in the legend

    Example:
        Create length distribution plots for a sample::

            fragments = Fragment.from_bam("sample.bam", "variants.vcf")
            fragment_length_plot(fragments, "plots/", "patient_001")
            # Creates: plots/patient_001_frags_len_kde.png

    See Also:
        fragment_length_gmm: For quantitative modeling of length distributions
    """
    logger: logging.Logger = get_logger()
    logger.info(
        "creating fragment length plots in ``{}``".format(out_dir)
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
        os.path.join(out_dir, "{}_frags_len_kde.png".format(name))

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
    """
    Fit Gaussian mixture models to fragment length distributions.

    Performs analysis of fragment length data by fitting Gaussian mixture
    models (GMMs) with a configurable number of components.

    The function automatically:

    * Fits GMM parameters using maximum likelihood estimation
    * Validates model convergence and goodness-of-fit
    * Generates diagnostic plots showing observed vs fitted distributions
    * Exports results as JSON, including parameters and statistics

    Args:
        fragments: Fragment list to analyze. All fragments are used
            except those marked as bogus (quality-filtered).
        config_filepath: Path to JSON configuration file specifying GMM
            parameters. Must contain number of components and initial values.
        out_dir: Output directory for results. Multiple files will be created
            including plots and parameter files.
        name: Base name for output files. Used to generate unique filenames
            for plots and parameter exports.

    Returns:
        None. Results are saved to multiple output files in the specified
            directory.

    Raises:
        SystemExit: If GMM fitting fails due to non-convergence or numerical
            issues. This indicates problematic data or inappropriate model
            specification.

    Note:
        * Requires a valid GMM configuration file (see ``configs/`` directory)
        * Generates both visual diagnostics and quantitative statistics
        * May not converge for poor data
        * Results include goodness-of-fit metrics for model validation

    Example:
        Fit a 3-component GMM to fragment lengths::

            fragments = Fragment.from_bam("sample.bam")
            fragment_length_gmm(
                fragments,
                "configs/gmm_3.json",
                "analysis/",
                "sample_001"
            )
            # Creates: analysis/sample_001_gmm_frags_len.json
            #          analysis/sample_001_gmm_plot.png

    See Also:
        * :func:`pyfraglib.math.fit_gmm` - Core GMM fitting algorithm
        * :func:`write_gmm_params` - Parameter export functionality
        * :doc:`../examples/index` - Complete analysis examples
    """
    logger: logging.Logger = get_logger()
    logger.info(
        "fitting GMM based on config file ``{}``, writing results to ``{}``"
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
    """
    Export Gaussian mixture model parameters and statistics to JSON.

    Serializes comprehensive GMM fitting results to a structured JSON file
    for downstream analysis, quality control, and reproducibility. The output
    includes fitted parameters, optimization metrics, and goodness-of-fit
    statistics.

    Args:
        n: Number of Gaussian components in the fitted model.
        params: Flattened array of fitted parameters in the order:
            [means, standard_deviations, mixing_weights]. Each component
            contributes one value to each parameter type.
        obj_val: Final objective function value from the optimization.
            Lower values indicate better fits for maximum likelihood
            estimation.
        conv: Boolean indicating whether the optimization algorithm converged
            successfully. False values suggest problematic fits.
        goodness_of_fit: Dictionary containing statistical measures of model
            quality, including Kolmogorov-Smirnov tests and Jensen-Shannon
            divergence.
        out_dir: Output directory path where the JSON file will be saved.
        name: Base filename. The output will be saved as
            ``{name}_gmm_frags_len.json``.

    Returns:
        None. Results are written to a JSON file in the specified directory.

    Note:
        * Parameters are automatically parsed and organized by component
        * Includes both optimization and statistical validation metrics
        * Mostly used as a utility function

    Example:
        Typical usage within GMM fitting workflow::

            # After fitting GMM
            write_gmm_params(
                n=3,
                params=[167, 120, 200, 25, 15, 30, 0.6, 0.3, 0.1],
                obj_val=-12345.67,
                conv=True,
                goodness_of_fit={"ks_statistic": 0.05, "ks_p_value": 0.8},
                out_dir="results/",
                name="sample"
            )
            # Creates: results/sample_gmm_frags_len.json

    File Format:
        The output JSON contains::

            {
                "number_of_gaussians": 3,
                "objective_value": -12345.67,
                "converged": true,
                "estimated_means": [167.0, 120.0, 200.0],
                "estimated_stds": [25.0, 15.0, 30.0],
                "estimated_pis": [0.6, 0.3, 0.1],
                "ks_statistic": 0.05,
                "ks_p_value": 0.8,
                ...
            }
    """
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
