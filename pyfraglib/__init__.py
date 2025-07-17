"""
``pyfraglib``: A Python Library for Cell-Free DNA Fragmentomics Analysis
========================================================================

``pyfraglib`` is a comprehensive Python library for analyzing high-throughput
sequencing data of cell-free DNA (cfDNA), with a focus on fragmentomics
features. It provides tools to extract, analyze, and visualize fragment
characteristics from BAM files, enabling researchers to study cfDNA biology and
develop liquid biopsy applications.

Key Features
------------
- **Fragment Extraction**: Extract cfDNA fragments from BAM files with quality
                           filtering and VCF-based mutation annotation
- **Batch Processing**: Efficient processing of large cohorts with *Nextflow*
- **Fragmentomics Analysis**: Calculate fragment lengths distributions, end
                              motifs, diversity and protection scores, ...
- **Advanced Statistical Analysis**: Gaussian mixture model fitting to fragment
                                     length distributions, non-negative matrix
                                     factorization, differential end motif
                                     abundance, ...
- **Visualization**: Generate publication-ready plots for fragment
                     characteristics (length plots, end motif plots, ...)
- **Simulation**: Generate synthetic cfDNA data for method validation
- **Additional tools**: Use CLI tools to subset BAM files or retrieve
                        transcription start side annotations

Core Classes
------------
- :class:`Fragment`: Represents a single cfDNA fragment with genomic
                     coordinates and properties
- :class:`FragmentList`: Collection of fragments from a single sample with
                         analysis methods
- :class:`FragmentCollection`: Multi-sample fragment collection for batch
                               processing
- :class:`FragFile`: Efficient reader/writer for serialized fragment data
                     (.frag files)
- :class:`FragmentSimulator`: Generate synthetic cfDNA fragments for method
                              validation
- :class:`TissueMixtureSimulator`: Simulate tissue mixtures for *in-silico*
                                   liquid biopsy experiments

Main Analysis Modules
---------------------
- :mod:`pyfraglib.fragment`: Core fragment classes and BAM file processing
- :mod:`pyfraglib.scores`: Fragmentomics scoring algorithms
- :mod:`pyfraglib.lengths`: Fragment length analysis and Gaussian mixture
                            modeling
- :mod:`pyfraglib.stats`: Overview statistics and visualization functions
- :mod:`pyfraglib.math`: Mathematical utilities for GMM fitting and statistics
- :mod:`pyfraglib.simulator`: Synthetic cfDNA data generation

Example Usage
-------------
Basic fragment analysis workflow:

.. code-block:: python

    from pyfraglib import Fragment, FragmentList
    from pyfraglib.lengths import fragment_length_plot
    from pyfraglib.scores import windowed_protection_score

    # Extract fragments from BAM file
    fragments = Fragment.from_bam("sample.bam", "variants.vcf")

    # Generate length distribution plot
    fragment_length_plot(fragments, "output/", "sample")

    # Calculate windowed protection score
    wps_scores = windowed_protection_score(fragments, "regions.bed")

    # Save fragments to file
    fragments.to_frag_file("sample", "output/")

Installation
------------
Install using pip:

.. code-block:: bash

    pip install pyfraglib

Or from source:

.. code-block:: bash

    git clone https://bitbucket.org/schwarzlab/pyfraglib.git
    cd pyfraglib
    pip install .

Command Line Interface
----------------------
``pyfraglib`` provides a command-line interface through the ``pyfrag.py``
script:

.. code-block:: bash

    # Display version information
    pyfrag.py version

    # Extract fragments from BAM files
    pyfrag.py -o output/ extract --bam-file sample.bam

    # Analyze fragment lengths
    pyfrag.py -o output/ lengths --frag-file sample.frag --config-file gmm.json

    # Calculate fragmentomics scores
    pyfrag.py -o output/ scores --frag-file sample.frag --bed-file regions.bed

    # Generate synthetic data
    pyfrag.py -o simulation/ simulate --config simulation.json

Requirements
------------
The following packages are required for ``pyfraglib`` to run:

- Python 3.10+
- NumPy
- SciPy
- Matplotlib
- Pandas
- pysam
- intervaltree
- tqdm

See ``setup.py`` for additional build dependencies.

Version Information
-------------------
Current version: Check ``pyfraglib.__version__``
Repository: https://bitbucket.org/schwarzlab/pyfraglib/

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

import logging

from pyfraglib.core import get_logger, fail, PyfragManager
from pyfraglib.math import fit_gmm, plot_gmm
from pyfraglib.fragment import Fragment, FragmentList, FragmentCollection
from pyfraglib.fragfile import FragFile
from pyfraglib.stats import fragments_per_chromosome_barplot
from pyfraglib.lengths import fragment_length_plot
from pyfraglib.simulator import FragmentSimulator, TissueMixtureSimulator, \
                                NucleaseProfile

__version__ = "0.5.1"
__repo_url__ = "https://bitbucket.org/schwarzlab/pyfraglib/"
__all__ = [
    "Fragment", "FragmentList", "FragmentCollection", "fail",
    "get_logger", "FragFile", "fragments_per_chromosome_barplot",
    "fragment_length_plot", "fit_gmm", "plot_gmm",
    "FragmentSimulator", "TissueMixtureSimulator", "NucleaseProfile"
]


class FixedLogLevelFilter(logging.Filter):
    """
    This filter changes the level of every record to ``level`` and only logs
    the record if that new level is more severe than the log level of the
    logger that emitted the record.
    """
    def __init__(self, level: int):
        super().__init__()
        self.level: int = level

    def filter(self, record: logging.LogRecord) -> bool:
        this_logger_name: str = record.name
        this_logger: logging.Logger = logging.getLogger(this_logger_name)

        record.levelno = self.level
        record.levelname = logging.getLevelName(self.level)
        return record.levelno >= this_logger.level


logging.basicConfig(
    level=logging.NOTSET,
    format='[%(asctime)s %(levelname)-8s %(name)-9s %(process)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logging.captureWarnings(True)

# @NOTE(ds): The ``matplotlib`` logging messages are just annoying. We turn
# them off completely.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
logging.getLogger("matplotlib.pyplot").setLevel(logging.ERROR)
logging.getLogger("PIL.PngImagePlugin").setLevel(logging.ERROR)

# @NOTE(ds): Not sure if the following log messages could be interesting for
# debugging.
logging.getLogger("asyncio").setLevel(logging.DEBUG)

# @NOTE(ds): ``py.warnings`` logs everything as a WARNING. We don't want that
# so we apply a filter to log everything as DEBUG.
py_warnings_logger: logging.Logger = logging.getLogger("py.warnings")
log_filter: FixedLogLevelFilter = FixedLogLevelFilter(logging.DEBUG)
py_warnings_logger.addFilter(log_filter)

PyfragManager.register(  # type: ignore
    "FragmentCollection", FragmentCollection
)
