Stats Command
=============

The ``stats`` command generates comprehensive statistics and visualizations from fragment files.

Syntax
------

.. code-block:: bash

   pyfrag.py stats [OPTIONS]

Options
-------

.. code-block:: bash

   --frag-file PATH       Single fragment file to analyze
   --frag-dir PATH        Directory containing fragment files
   --out-dir PATH         Output directory for plots and statistics (required)

Examples
--------

Analyze Single Fragment File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py stats --frag-file sample.frag --out-dir statistics/

Analyze Multiple Fragment Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py stats --frag-dir fragments/ --out-dir statistics/

Output
------

The stats command generates several output files:

Plots
~~~~~

- ``*_fragment_stats.png`` - Fragment length distribution histogram
- ``*_chromosome_distribution.png`` - Fragments per chromosome bar plot
- ``*_end_motifs_5p.png`` - 5' end motif frequency bar plot
- ``*_end_motifs_3p.png`` - 3' end motif frequency bar plot

Statistics Files
~~~~~~~~~~~~~~~~

- ``*_summary.txt`` - Basic fragment statistics
- ``*_motif_counts.csv`` - Detailed end motif counts
- ``*_length_stats.csv`` - Fragment length statistics by chromosome

Generated Statistics
--------------------

Basic Fragment Statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~

- Total fragment count
- Bogus fragment count and percentage
- Mutated fragment count and percentage (if VCF was used)
- Mean and median fragment length
- Fragment length standard deviation

Chromosome Distribution
~~~~~~~~~~~~~~~~~~~~~~~

- Fragment count per chromosome
- Percentage distribution across chromosomes
- Length statistics per chromosome

End Motif Analysis
~~~~~~~~~~~~~~~~~~

- 3-mer and 4-mer frequency analysis
- Shannon entropy and Simpson diversity index
- Most common motif sequences
- Motif diversity comparison between 5' and 3' ends

Requirements
------------

- Fragment files must be valid .frag files created by the extract command
- Sufficient disk space for output plots and statistics
- matplotlib for plot generation

Performance Notes
-----------------

- Statistics generation is memory-efficient
- Processing time scales with fragment count
- Plot generation may take longer for large datasets
- Multiple files are processed sequentially

Troubleshooting
---------------

**"Fragment file corrupted" error**
   Re-extract fragments from the original BAM file

**Memory errors with large files**
   Process files individually rather than in batches

**Missing plots**
   Check that matplotlib is properly installed and configured