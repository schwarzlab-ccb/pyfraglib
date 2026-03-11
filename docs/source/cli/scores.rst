Scores Command
==============

The ``scores`` command calculates fragmentomics scores like the Windowed Protection Score (WPS) and motif diversity metrics.

Syntax
------

.. code-block:: bash

   pyfrag.py scores -o <OUT_DIR> [OPTIONS]

Options
-------

.. code-block:: bash

   --frag-file PATH       Single fragment file to analyze
   --frag-dir PATH        Directory containing fragment files
   --bed-file PATH        BED file defining genomic regions (required)

Examples
--------

.. code-block:: bash

   pyfrag.py --out-dir scores/ scores --frag-file sample.frag --bed-file regions.bed


BED File Format
---------------

The BED file should contain genomic regions of interest for which the windowed protection score is then calculated:

.. code-block:: text

   chr1    1000000    1001000    region1
   chr1    2000000    2001000    region2
   chr2    500000     501000     region3

Required columns:
- Column 1: Chromosome name
- Column 2: Start position (0-based)
- Column 3: End position (0-based)
- Column 4: Region name (optional)

Output
------

As output, this subcommand saves *CSV* files containing the windowed protection score per position as well as global motif diversity metrics to disk. A genome-wide WPS plot is saved, too.

Calculated Metrics
------------------

Windowed Protection Score (WPS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For every genomic position :math:`i`, the WPS measures the depletion of fragment ends within a window surrounding that position:

.. math::

   \text{WPS}(i) = N_{span}(i) - N_{end}(i)



Large values indicate protection (fewer fragment ends), whereas smaller values indicate chromatin accessibility (more fragment ends). Because this metric is not independent of local fragment depth, the output *CSV* file contains positional fragment depth, the number of ending, and the number of spanning fragments as additional fields. This way, users can calculate fractions instead of differences if that makes more sense for their datasets.

Motif Diversity
~~~~~~~~~~~~~~~

Given an end motif distribution :math:`X = (x_1, ..., x_n)` where :math:`x_i` is the observed count for end motif :math:`i`, we calculate Shannon entropy :math:`H(X)` as

.. math::

   H(X) = -\sum_{i=1}^{n} p(x_i) \log_2 p(x_i)

and Simpson index :math:`D(X)` as

.. math::

    D(X) = \sum_{i=1}^{n} x_i^2

.
