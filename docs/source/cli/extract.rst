Extract Command
===============

The ``extract`` command processes BAM files to extract cfDNA fragments and optionally annotate them with mutation information from VCF files.

Syntax
------

.. code-block:: bash

   pyfrag.py -o <OUT_DIR> extract [OPTIONS]

The ``-o / --out-dir`` must come *before* the subcommand. Please be aware that all other subcommands follow this pattern, too.

Options
-------

.. code-block:: bash

   --bam-file PATH        Single BAM file to process
   --bam-dir PATH         Directory containing BAM files to process
   --with-vcf             Search for VCF file for mutation annotation (optional)
   --nanopore             Process BAM as single-ended Nanopore data

If the ``--with-vcf`` flag is set, ``pyfraglib`` will search for one *VCF* file per *BAM*. The *VCF* must be named identically, except for the file extension which must be ``.vcf`` instead of ``.bam``. *BAM* files must be indexed.

Examples
--------

.. code-block:: bash

   pyfrag.py --out-dir fragments/ extract --bam-file sample.bam --with-vcf


Output
------

The extract command creates ``.frag`` files in the specified output directory. These files contain compressed, serialized fragment data that can be used by other pyfraglib commands. For an input file named ``sample.bam``, an output file called ``sample.frag`` will be created in ``<OUT_DIR>``.

Quality Filters
---------------

The extraction process applies several quality filters that are currently hard-coded as source code constants:

- **MAPQ ≥ 20**: Minimum mapping quality score
- **Insert size ≤ 900bp**: Maximum insert size for paired-end reads
- **Valid chromosomes**: Only standard chromosomes (1-22, X, Y, M)
- **No 'N' bases**: Filters out fragments with ambiguous nucleotides
- **Proper pairing**: For paired-end data, both reads must be properly paired

General Tips
------------

If you want to process a large set of samples, we recommend against using the ``--bam-dir`` flag. Even though the latter works on *BAM* files in parallel and writes fragment data to disk immediately, it still requires potentially large amounts of memory. We much rather recommend adopting the Nextflow pipeline in our ``tools/`` directory for working with many samples.
