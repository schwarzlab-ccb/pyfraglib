Utility Scripts
===============

pyfraglib includes several utility scripts for specialized tasks. While they were created for solving specific problems during our research, they might still be useful to others, so we briefly describe them here.

extract_mutated_reads.py
========================

Extracts reads that support specific variants from BAM files. This script can be used to split a *BAM* file into separate files containing only the mutated reads vs. only the unmutated reads based on a *VCF* file. Mutational status assignment in ``pyfraglib`` was extensively tested using *BAM* files prepared in this way.

Syntax
------

.. code-block:: bash

   extract_mutated_reads.py [OPTIONS]

Options
-------

.. code-block:: bash

   --bam PATH        Input BAM file
   --vcf PATH        Input VCF file with variants
   --output PATH     Output BAM file for mutated reads
   --keep            Keep unmutated reads in separate BAM file


download_tss_annos.py
=====================

Downloads transcription start site (TSS) annotations from Ensembl, generating a *BED* file that can be used with the ``pyfraglib scores`` subcommand. The genes for which to download annotations must be provided via a *TXT* file with one gene per line.

Syntax
------

.. code-block:: bash

   download_tss_annos.py [OPTIONS]

Options
-------

.. code-block:: bash

   --ref VERSION          Genome version (hg19, hg38)
   --outfile PATH         Output directory for annotations
   --gene-list PATH       File with gene symbols to download
   --verbose              Enable verbose logging (e.g. for trouble-shooting)


txt_to_vcf.py
=============

Converts custom variant format to standard VCF format. This script is probably not useful for anyone but researchers from our group.

Syntax
------

.. code-block:: bash

   txt_to_vcf.py [OPTIONS]

Options
-------

.. code-block:: bash

   --infile PATH      Input file in custom format
   --outfile PATH     Output VCF file
   --ref-genome PATH  hg19 or hg38

