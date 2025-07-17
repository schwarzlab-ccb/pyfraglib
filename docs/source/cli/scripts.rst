Utility Scripts
===============

pyfraglib includes several utility scripts for specialized tasks.

extract_mutated_reads.py
========================

Extracts reads that support specific variants from BAM files.

Syntax
------

.. code-block:: bash

   python scripts/extract_mutated_reads.py [OPTIONS]

Options
-------

.. code-block:: bash

   --bam-file PATH        Input BAM file
   --vcf-file PATH        Input VCF file with variants
   --output-bam PATH      Output BAM file for mutated reads
   --keep                 Keep unmutated reads in separate BAM file

Examples
--------

Extract Mutated Reads Only
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python scripts/extract_mutated_reads.py --bam-file sample.bam --vcf-file variants.vcf --output-bam mutated_reads.bam

Keep Both Mutated and Unmutated Reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python scripts/extract_mutated_reads.py --bam-file sample.bam --vcf-file variants.vcf --output-bam mutated_reads.bam --keep

Output
------

- ``mutated_reads.bam`` - Reads supporting variant alleles
- ``unmutated_reads.bam`` - Reads with reference alleles (if --keep used)
- Index files (.bai) for both output files

download_tss_annos.py
=====================

Downloads transcription start site (TSS) annotations from Ensembl.

Syntax
------

.. code-block:: bash

   python scripts/download_tss_annos.py [OPTIONS]

Options
-------

.. code-block:: bash

   --genome VERSION       Genome version (hg19, hg38)
   --output-dir PATH      Output directory for annotations
   --gene-list PATH       File with gene symbols to download (optional)
   --upstream INT         Upstream region size (default: 1000)
   --downstream INT       Downstream region size (default: 1000)

Examples
--------

Download All TSS Annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python scripts/download_tss_annos.py --genome hg38 --output-dir annotations/

Download Specific Genes
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python scripts/download_tss_annos.py --genome hg19 --gene-list gene_list.txt --output-dir annotations/

Custom Region Sizes
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python scripts/download_tss_annos.py --genome hg38 --upstream 2000 --downstream 500 --output-dir annotations/

Output
------

- ``tss_annotations_{genome}.bed`` - BED file with TSS regions
- ``gene_annotations_{genome}.tsv`` - Detailed gene information
- ``download_log.txt`` - Download log and statistics

txt_to_vcf.py
=============

Converts custom variant format to standard VCF format.

Syntax
------

.. code-block:: bash

   python scripts/txt_to_vcf.py [OPTIONS]

Options
-------

.. code-block:: bash

   --input-file PATH      Input file in custom format
   --output-file PATH     Output VCF file
   --reference PATH       Reference FASTA file
   --sample-name NAME     Sample name for VCF header

Examples
--------

Basic Conversion
~~~~~~~~~~~~~~~~

.. code-block:: bash

   python scripts/txt_to_vcf.py --input-file variants.txt --output-file variants.vcf --reference reference.fasta --sample-name SAMPLE1

Input Format
~~~~~~~~~~~~

The input file should contain tab-separated values:

.. code-block:: text

   chr1    1000000    A    T    0.3
   chr1    2000000    G    C    0.15
   chr2    500000     C    A    0.25

Columns:
1. Chromosome
2. Position (1-based)
3. Reference allele
4. Alternative allele
5. Allele frequency

Output
------

- Standard VCF file with proper headers
- Validation log with conversion statistics
- Error log for problematic variants

Common Workflows
================

Complete Fragment Analysis Pipeline
-----------------------------------

.. code-block:: bash

   # 1. Extract fragments
   pyfrag.py extract --bam-file sample.bam --vcf-file variants.vcf --out-dir fragments/
   
   # 2. Generate statistics
   pyfrag.py stats --frag-file fragments/sample.frag --out-dir stats/
   
   # 3. Analyze lengths
   pyfrag.py lengths --frag-file fragments/sample.frag --config-file configs/gmm_3.json --out-dir lengths/
   
   # 4. Calculate scores
   pyfrag.py scores --frag-file fragments/sample.frag --bed-file tss_regions.bed --out-dir scores/

Mutation Analysis Workflow
---------------------------

.. code-block:: bash

   # 1. Extract mutated reads
   python scripts/extract_mutated_reads.py --bam-file sample.bam --vcf-file variants.vcf --output-bam mutated.bam --keep
   
   # 2. Process both BAM files
   pyfrag.py extract --bam-file mutated.bam --out-dir fragments/
   pyfrag.py extract --bam-file unmutated.bam --out-dir fragments/
   
   # 3. Compare fragment characteristics
   pyfrag.py stats --frag-file fragments/mutated.frag --out-dir stats/mutated/
   pyfrag.py stats --frag-file fragments/unmutated.frag --out-dir stats/unmutated/

Simulation Validation Workflow
-------------------------------

.. code-block:: bash

   # 1. Generate synthetic data
   pyfrag.py simulate --config configs/simulation_example.json --out-dir synthetic/
   
   # 2. Analyze synthetic fragments
   pyfrag.py stats --frag-file synthetic/synthetic_sample.frag --out-dir synthetic/stats/
   
   # 3. Compare with real data
   pyfrag.py stats --frag-file fragments/real_sample.frag --out-dir real/stats/

Troubleshooting
===============

Common Issues
-------------

**Script not found errors**
   Ensure you're running scripts from the pyfraglib root directory

**Permission errors**
   Check file permissions and output directory write access

**Memory errors**
   Reduce batch sizes or process files individually

**Format errors**
   Validate input file formats before processing

Getting Help
------------

For help with any script:

.. code-block:: bash

   python scripts/script_name.py --help

Or consult the API documentation for programmatic usage.