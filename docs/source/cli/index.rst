Command Line Interface
======================

pyfraglib provides a comprehensive command-line interface through the ``pyfrag.py`` script. This section covers all available commands and their options.

.. toctree::
   :maxdepth: 2

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
   --version      Show program version

Getting Help
------------

For help with any subcommand:

.. code-block:: bash

   pyfrag.py <subcommand> --help

Example: ``pyfrag.py extract --help``

Basic Usage Examples
--------------------

Extract fragments from a BAM file:

.. code-block:: bash

   pyfrag.py extract --bam-file sample.bam --out-dir fragments/

Generate statistics:

.. code-block:: bash

   pyfrag.py stats --frag-file sample.frag --out-dir statistics/

Calculate WPS scores:

.. code-block:: bash

   pyfrag.py scores --frag-file sample.frag --bed-file regions.bed --out-dir scores/

Simulate synthetic data:

.. code-block:: bash

   pyfrag.py simulate --config simulation_config.json --out-dir synthetic/

Additional Scripts
------------------

pyfraglib also includes several utility scripts:

* ``extract_mutated_reads.py`` - Extract reads supporting variants
* ``download_tss_annos.py`` - Download TSS annotations
* ``txt_to_vcf.py`` - Convert custom format to VCF