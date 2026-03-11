Lengths Command
===============

The ``lengths`` command plots fragment length distributions and fits a configurable Gaussian Mixture Model (GMM) to the distributions. Model parameters and visualizations are saved to disk.

Syntax
------

.. code-block:: bash

   pyfrag.py lengths -o <OUT_DIR> [OPTIONS]

Options
-------

.. code-block:: bash

   --frag-file PATH       Single fragment file to analyze
   --frag-dir PATH        Directory containing fragment files
   --config-file PATH     GMM configuration file (JSON format, required)

Examples
--------

.. code-block:: bash

   pyfrag.py --out-dir lengths/ lengths --frag-file sample.frag --config-file configs/gmm_3.json

Configuration File Format
--------------------------

GMM configuration files use *JSON* format with the following fields:

.. code-block:: json

    {
        "number_of_gaussians": 3,
        "subsample_percentage": 0.10,
        "means_lower_bounds": [50, 250, 450],
        "means_upper_bounds": [250, 450, 650],
        "std_lower_bounds": [1, 1, 1],
        "std_upper_bounds": [30, 60, 90],
        "initial_means": [160, 320, 480],
        "initial_pis": [0.70, 0.20, 0.10]
    }

The parameters determine the number of Gaussians fit, the lower and upper bounds on the means and standard deviations, and the initial parameters for means, standard deviations, and mixture fractions (:math:`\pi`). The `subsample_percentage` parameter can be used if the `.frag` file contains a large number of fragments which sometimes makes the GMM fitting prohibitively slow. The model indicates such situations with an informative error message.

Output
------

The output saved to ``<OUT_DIR>`` contains a kernel density estimate plot of the fragment length distribution, stratified into mutated vs. wildtype fragments, and a plot of the GMM fit. Model parameters are saved to a *JSON* file along with goodness-of-fit parameters.

