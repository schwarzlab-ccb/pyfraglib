Lengths Command
===============

The ``lengths`` command performs Gaussian Mixture Model (GMM) analysis on fragment length distributions.

Syntax
------

.. code-block:: bash

   pyfrag.py lengths [OPTIONS]

Options
-------

.. code-block:: bash

   --frag-file PATH       Single fragment file to analyze
   --frag-dir PATH        Directory containing fragment files
   --config-file PATH     GMM configuration file (JSON format)
   --out-dir PATH         Output directory for plots and results (required)

Examples
--------

Analyze Single Fragment File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py lengths --frag-file sample.frag --config-file configs/gmm_3.json --out-dir lengths/

Analyze Multiple Fragment Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py lengths --frag-dir fragments/ --config-file configs/gmm_2.json --out-dir lengths/

Configuration File Format
--------------------------

GMM configuration files use JSON format:

.. code-block:: json

   {
       "n_components": 3,
       "means_init": [167, 320, 450],
       "covariance_type": "full",
       "max_iter": 100,
       "n_init": 10,
       "random_state": 42
   }

Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

- ``n_components``: Number of Gaussian components to fit
- ``means_init``: Initial means for components (optional)
- ``covariance_type``: Type of covariance matrix ("full", "tied", "diag", "spherical")
- ``max_iter``: Maximum number of iterations for fitting
- ``n_init``: Number of random initializations
- ``random_state``: Random seed for reproducibility

Output
------

Plots
~~~~~

- ``*_length_distribution.png`` - Histogram with GMM overlay
- ``*_gmm_components.png`` - Individual component plots
- ``*_residuals.png`` - Residual analysis plots

Results Files
~~~~~~~~~~~~~

- ``*_gmm_results.json`` - Fitted parameters and statistics
- ``*_components.csv`` - Component means, weights, and covariances
- ``*_goodness_of_fit.csv`` - AIC, BIC, and log-likelihood values

Interpretation
--------------

Component Analysis
~~~~~~~~~~~~~~~~~~

Typical cfDNA shows multiple fragment size components:

- **Nucleosomal fragments** (~167bp): Mono-nucleosomal DNA
- **Di-nucleosomal fragments** (~320bp): Di-nucleosomal DNA  
- **Longer fragments** (>400bp): Higher-order chromatin structures

Goodness of Fit
~~~~~~~~~~~~~~~

- **AIC** (Akaike Information Criterion): Lower values indicate better fit
- **BIC** (Bayesian Information Criterion): Penalizes model complexity
- **Log-likelihood**: Higher values indicate better fit

Requirements
------------

- Fragment files must be valid .frag files
- Configuration file must be valid JSON
- scipy for GMM fitting
- matplotlib for plotting

Performance Notes
-----------------

- GMM fitting is computationally intensive
- Memory usage scales with fragment count
- Multiple components increase fitting time
- Convergence may be slow for complex distributions

Troubleshooting
---------------

**"GMM failed to converge" error**
   Try different initialization parameters or increase max_iter

**Poor fit quality**
   Adjust the number of components or initialization means

**Memory errors**
   Reduce the number of components or process smaller datasets