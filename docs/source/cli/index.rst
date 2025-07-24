Command Line Interface
======================

pyfraglib provides a command-line interface through the ``pyfrag.py`` script. This section covers all available commands and most of their options.

.. toctree::
   :maxdepth: 1

   extract
   stats
   lengths
   scores
   simulate
   scripts

Overview
--------

The main CLI script is ``pyfrag.py``, which provides several subcommands:

.. code-block:: bash

   pyfrag.py <subcommand> [options]

Available subcommands:

* ``extract`` - Extract fragments from BAM files
* ``stats`` - Generate fragment statistics
* ``lengths`` - Analyze fragment length distributions
* ``scores`` - Calculate fragmentomics scores
* ``simulate`` - Generate synthetic cfDNA data
* ``version`` - Show version information

Global Options
--------------

All subcommands support these global options:

.. code-block:: bash

   -h, --help     Show help message and exit
   -v, --verbose  Enable verbose logging
   -o, --out-dir  Set the output directory for the subcommand

Importantly, for help with any subcommand:

.. code-block:: bash

   pyfrag.py <subcommand> --help


Additional Scripts
------------------

pyfraglib also includes several utility scripts:

* ``extract_mutated_reads.py`` - Extract reads supporting variants
* ``download_tss_annos.py`` - Download TSS annotations
* ``txt_to_vcf.py`` - Convert custom format to VCF
