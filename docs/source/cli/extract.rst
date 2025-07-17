Extract Command
===============

The ``extract`` command processes BAM files to extract cfDNA fragments and optionally annotate them with mutation information from VCF files.

Syntax
------

.. code-block:: bash

   pyfrag.py extract [OPTIONS]

Options
-------

Input Options
~~~~~~~~~~~~~

.. code-block:: bash

   --bam-file PATH        Single BAM file to process
   --bam-dir PATH         Directory containing BAM files
   --vcf-file PATH        VCF file for mutation annotation (optional)
   --vcf-dir PATH         Directory containing VCF files (optional)

Output Options
~~~~~~~~~~~~~~

.. code-block:: bash

   --out-dir PATH         Output directory for .frag files (required)

Processing Options
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   --nanopore             Process as single-ended Nanopore data
   --parallel             Enable parallel processing for multiple files

Examples
--------

Process Single BAM File
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py extract --bam-file sample.bam --out-dir fragments/

Process with Mutation Annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py extract --bam-file sample.bam --vcf-file variants.vcf --out-dir fragments/

Process Multiple BAM Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py extract --bam-dir bam_files/ --out-dir fragments/

Process Nanopore Data
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py extract --bam-file nanopore.bam --nanopore --out-dir fragments/

Output
------

The extract command creates ``.frag`` files in the specified output directory. These files contain compressed, serialized fragment data that can be used by other pyfraglib commands.

Output file naming:
- Input: ``sample.bam`` → Output: ``sample.frag``
- Input: ``data/patient1.bam`` → Output: ``patient1.frag``

Quality Filters
---------------

The extraction process applies several quality filters:

- **MAPQ ≥ 20**: Minimum mapping quality score
- **Insert size ≤ 900bp**: Maximum insert size for paired-end reads
- **Valid chromosomes**: Only standard chromosomes (1-22, X, Y, M)
- **No 'N' bases**: Filters out fragments with ambiguous nucleotides
- **Proper pairing**: For paired-end data, both reads must be properly paired

Requirements
------------

- BAM files must be indexed (``.bai`` files required)
- Sufficient disk space for output files
- VCF files must match the reference genome used for BAM alignment

Performance Tips
----------------

- Use SSD storage for better I/O performance
- Process multiple files in parallel when possible
- Ensure adequate RAM (8GB+ recommended for large BAM files)
- Pre-sort and index BAM files for optimal performance

Troubleshooting
---------------

**"BAM file not indexed" error**
   Create an index with: ``samtools index input.bam``

**Memory errors**
   Reduce the number of files processed simultaneously or increase available RAM

**No fragments extracted**
   Check that BAM file contains properly aligned reads with valid chromosome names