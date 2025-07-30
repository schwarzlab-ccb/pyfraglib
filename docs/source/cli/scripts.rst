Utility Scripts
===============

pyfraglib includes several utility scripts for specialized tasks. While they were created for solving specific problems during our research, they might still be useful to others, so we briefly describe them here.

extract_mutated_reads.py
------------------------

Extracts reads that support specific variants from BAM files. This script can be used to split a *BAM* file into separate files containing only the mutated reads vs. only the unmutated reads based on a *VCF* file. Mutational status assignment in ``pyfraglib`` was extensively tested using *BAM* files prepared in this way.

Syntax
~~~~~~

.. code-block:: bash

   extract_mutated_reads.py [OPTIONS]

Options
~~~~~~~

.. code-block:: bash

   --bam PATH        Input BAM file
   --vcf PATH        Input VCF file with variants
   --output PATH     Output BAM file for mutated reads
   --keep            Keep unmutated reads in separate BAM file


download_tss_annos.py
---------------------

Downloads transcription start site (TSS) annotations from Ensembl, generating a *BED* file that can be used with the ``pyfraglib scores`` subcommand. The genes for which to download annotations must be provided via a *TXT* file with one gene per line.

Syntax
~~~~~~

.. code-block:: bash

   download_tss_annos.py [OPTIONS]

Options
~~~~~~~

.. code-block:: bash

   --ref VERSION          Genome version (hg19, hg38)
   --outfile PATH         Output directory for annotations
   --gene-list PATH       File with gene symbols to download
   --verbose              Enable verbose logging (e.g. for trouble-shooting)


txt_to_vcf.py
-------------

Converts custom variant format to standard VCF format. This script is probably not useful for anyone but researchers from our group.

Syntax
~~~~~~

.. code-block:: bash

   txt_to_vcf.py [OPTIONS]

Options
~~~~~~~

.. code-block:: bash

   --infile PATH      Input file in custom format
   --outfile PATH     Output VCF file
   --ref-genome PATH  hg19 or hg38


nmf_fragment_lengths.py
-----------------------

Performs non-negative matrix factorization (NMF) on fragment length distributions from cfDNA samples. This analysis identifies underlying signatures in fragment length patterns that can reveal tissue-specific or pathological fragmentation patterns.

NMF decomposes the fragment length matrix into basis components (signatures) and mixing coefficients (sample compositions). The script runs multiple random initializations to ensure stable results and normalizes mixing coefficients to represent true proportions.

Input Requirements
~~~~~~~~~~~~~~~~~~

* Text file listing fragment length CSV file paths (one per line)
* CSV files from ``pyfrag stats`` with columns: ``fragment_length``, ``count``

Output Files
~~~~~~~~~~~~

* ``nmf_mixing_coefficients.csv`` - Sample composition weights (sum to 1.0 per sample)
* ``nmf_signatures.csv`` - Component signatures across fragment lengths
* ``nmf_signatures.png`` - Component signature line plots
* ``nmf_sample_composition.png`` - Sample composition heatmap
* ``nmf_component_weights.png`` - Component weights bar chart

Syntax
~~~~~~

.. code-block:: bash

   nmf_fragment_lengths.py [OPTIONS]

Options
~~~~~~~

.. code-block:: bash

   --file-list PATH           Text file with CSV file paths (one per line)
   --n-components INTEGER     Number of NMF components to extract (default: 3)
   --out-dir PATH            Output directory (default: nmf_output)
   --verbose                 Enable verbose logging


differential_end_motifs.py
--------------------------

Performs differential analysis of fragment end motif abundances between two groups of cfDNA samples. This statistical analysis identifies motifs that are significantly over-represented in one group versus another, which can reveal tissue-specific cleavage preferences or pathological changes in nuclease activity.

The analysis uses Wilcoxon rank-sum tests for non-parametric comparison, followed by Benjamini-Hochberg FDR correction to control for multiple testing across all motifs.

Input Requirements
~~~~~~~~~~~~~~~~~~

* JSON configuration file with ``group_a`` and ``group_b`` fields
* Each group contains a list of CSV file paths
* CSV files from ``pyfrag stats`` with columns: ``motif_5p``, ``count_5p``, ``motif_3p``, ``count_3p``

Output Files
~~~~~~~~~~~~

* ``differential_results.csv`` - Complete statistical results for all motifs
* ``differential_volcano_plot.png`` - Effect size vs significance visualization
* ``differential_effect_sizes.png`` - Distribution of effect sizes
* ``differential_top_motifs.png`` - Top significant motifs bar chart

Syntax
~~~~~~

.. code-block:: bash

   differential_end_motifs.py [OPTIONS]

Options
~~~~~~~

.. code-block:: bash

   --config PATH         JSON configuration file with group definitions
   --out-dir PATH        Output directory (default: differential_output)
   --verbose             Enable verbose logging

Configuration Format
~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

   {
     "group_a": [
       "healthy1_k4_end_motifs.csv",
       "healthy2_k4_end_motifs.csv",
       "healthy3_k4_end_motifs.csv"
     ],
     "group_b": [
       "cancer1_k4_end_motifs.csv",
       "cancer2_k4_end_motifs.csv",
       "cancer3_k4_end_motifs.csv"
     ]
   }

