Installation
============

``pyfraglib`` depends on several external Python packages. The recommended installation method below is ``setup.py`` and ``pyproject.toml`` based so dependencies should be installed for you. For ``pysam``, please refer to their documentation to learn about common issues (if you have any, that is!). Please also take a look at our ``README.md`` for the latest installation instructions. Installing ``pyfraglib`` from PyPI will be supported in the future.

Conda Installation (Recommended)
---------------------------------

1. Create the conda environment from the provided YAML file:

.. code-block:: bash

   conda env create -f pyfraglib.yml
   conda activate pyfraglib

2. Install pyfraglib into the environment:

.. code-block:: bash

   python -m pip install .

Development Installation
------------------------

For development work, use the development installation script:

.. code-block:: bash

   ./tools/dev_install.sh

This script will run type checking with mypy, perform linting with flake8, install the package, and run all tests.

Verification
------------

To verify your installation, run the test suite:

.. code-block:: bash

   python -m unittest discover -s tests -v

Some integration tests might be reported as skipped if the necessary data files are not available. We do not provide data files for our more complex tests for privacy reasons.

You can have ``pyfraglib`` print version information as follows:

.. code-block:: bash

   pyfrag.py version

