Scores Command
==============

The ``scores`` command calculates fragmentomics scores including Windowed Protection Score (WPS) and motif diversity metrics.

Syntax
------

.. code-block:: bash

   pyfrag.py scores [OPTIONS]

Options
-------

.. code-block:: bash

   --frag-file PATH       Single fragment file to analyze
   --frag-dir PATH        Directory containing fragment files
   --bed-file PATH        BED file defining genomic regions (required)
   --out-dir PATH         Output directory for scores and plots (required)
   --window-size INT      Window size for WPS calculation (default: 120)
   --step-size INT        Step size for sliding window (default: 5)

Examples
--------

Calculate WPS for Single Sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py scores --frag-file sample.frag --bed-file regions.bed --out-dir scores/

Calculate Scores for Multiple Samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py scores --frag-dir fragments/ --bed-file tss_regions.bed --out-dir scores/

Custom Window Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py scores --frag-file sample.frag --bed-file regions.bed --window-size 100 --step-size 10 --out-dir scores/

BED File Format
---------------

The BED file should contain genomic regions of interest:

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

Score Files
~~~~~~~~~~~

- ``*_wps_scores.csv`` - WPS scores for each region and position
- ``*_motif_diversity.csv`` - Motif diversity metrics per region
- ``*_region_summary.csv`` - Summary statistics per region

Plots
~~~~~

- ``*_wps_profile.png`` - WPS profile plots for each region
- ``*_motif_diversity.png`` - Motif diversity bar plots
- ``*_score_heatmap.png`` - Heatmap of scores across regions

Calculated Metrics
------------------

Windowed Protection Score (WPS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WPS measures the relative depletion of fragment ends within a window compared to flanking regions:

.. math::

   WPS = \frac{Coverage_{flanking} - Coverage_{window}}{Coverage_{flanking} + Coverage_{window}}

Properties:
- Range: -1 to +1
- Positive values indicate protection (fewer fragment ends)
- Negative values indicate accessibility (more fragment ends)

Motif Diversity
~~~~~~~~~~~~~~~

Calculated for 3-mer and 4-mer end motifs:

- **Shannon Entropy**: Measures motif diversity
- **Simpson Index**: Measures motif dominance
- **Evenness**: Measures motif distribution uniformity

Requirements
------------

- Fragment files must be valid .frag files
- BED file must be properly formatted
- Genomic regions should be reasonable size (1-10kb recommended)
- Sufficient memory for score calculations

Performance Notes
-----------------

- WPS calculation is computationally intensive
- Memory usage scales with region size and fragment count
- Parallel processing used for multiple regions
- Large regions may require significant processing time

Biological Interpretation
-------------------------

WPS Applications
~~~~~~~~~~~~~~~~

- **Transcription Start Sites (TSS)**: High WPS indicates active transcription
- **Nucleosome Positioning**: WPS patterns reveal nucleosome organization
- **Chromatin Accessibility**: Low WPS suggests open chromatin
- **Regulatory Elements**: WPS changes indicate regulatory activity

Motif Diversity Applications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Fragmentation Patterns**: Different tissues show distinct motif patterns
- **Nuclease Activity**: Motif diversity reflects nuclease preferences
- **Sample Quality**: Low diversity may indicate degradation

Troubleshooting
---------------

**"No fragments in region" error**
   Check that BED regions overlap with fragment data

**Memory errors**
   Use smaller regions or reduce the number of regions processed simultaneously

**Poor WPS signal**
   Ensure adequate fragment coverage in target regions