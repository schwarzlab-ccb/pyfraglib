Stats Command
=============

The ``stats`` command generates a basic set of statistics and visualizations from fragment files.

Syntax
------

.. code-block:: bash

   pyfrag.py -o <OUT_DIR> stats [OPTIONS]

Options
-------

.. code-block:: bash

   --frag-file PATH       Single fragment file to analyze
   --frag-dir PATH        Directory containing fragment files

Examples
--------

.. code-block:: bash

   pyfrag.py stats --frag-file sample.frag --out-dir statistics/

Output
------

The stats command generates plots and data files:

- ``*_frag_stats.json``: summary statistics regarding the extracted fragments like the total and mutated number of fragments
- ``*_mut_frags_per_chrom.png``: bar plot of mutated vs. wildtype fragments per chromosome
- ``**_Xk_Pp_frag_end_motifs.png``: count plots of *X*-mer 5' / 3' end motifs (*X* is currently not parameterized but set to a constant *4*)
