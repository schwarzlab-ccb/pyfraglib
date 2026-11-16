Extract Command
===============

The ``extract`` command processes BAM files to extract cfDNA fragments and optionally annotate them with mutation information from VCF files.

Syntax
------

.. code-block:: bash

   pyfrag.py extract -o <OUT_DIR> [OPTIONS]

Options
-------

.. code-block:: bash

   --bam-file PATH          Single BAM file to process
   --bam-dir PATH           Directory containing BAM files to process
   --with-vcf               Search for VCF file for mutation annotation (optional)
   --nanopore               Process BAM as single-ended (e.g. Nanopore) data
   --skip-supplementary     Ignore secondary and supplementary alignments
   --max-insert-size INT    Upper bound on fragment length in bp (default: 900)
   --min-mapq INT           Minimum mapping quality per read (default: 20)

If the ``--with-vcf`` flag is set, ``pyfraglib`` will search for one *VCF* file per *BAM*. The *VCF* must be named identically, except for the file extension which must be ``.vcf`` instead of ``.bam``. *BAM* files must be indexed.

Examples
--------

.. code-block:: bash

   # paired-end short reads with mutation annotation
   pyfrag.py extract --bam-file sample.bam --with-vcf --out-dir fragments/

   # single-end long reads (e.g. ONT)
   pyfrag.py extract --bam-file sample.bam --nanopore --skip-supplementary \
       --max-insert-size 5000 --out-dir fragments/

Single-End (Nanopore) Data
--------------------------

With ``--nanopore``, every read is treated as one fragment spanning its aligned reference coordinates. Soft-clipped bases (e.g. adapters and barcodes) are ignored.

A single-end read only observes the fragment end at which sequencing starts reliably; the end where the read stops is prone to base-calling artefacts. For end motifs, ``pyfraglib`` therefore uses only the end where sequencing started: 5' end motifs come from forward-strand reads and 3' end motifs from reverse-strand reads. As for paired-end data, both are reported in reference orientation. The other end of each read is stored as unknown (``NNNN``) and is not counted in end motif analyses. Since both strands are sequenced at similar rates, 5' and 3' motif distributions each remain unbiased, but are based on roughly half of the fragments. Fragment lengths and coordinates use both ends.

Long-read aligners frequently report split (supplementary) and multi-mapping (secondary) alignments. By default, these are counted as additional fragments, but setting ``--skip-supplementary`` is recommended for long-read data so that each molecule yields at most one fragment. Consider raising ``--max-insert-size`` as well, since long reads span fragments far beyond the paired-end default of 900 bp.


Output
------

The extract command creates ``.frag`` files in the specified output directory. These files contain compressed, serialized fragment data that can be used by other pyfraglib commands. For an input file named ``sample.bam``, an output file called ``sample.frag`` will be created in ``<OUT_DIR>``.

Quality Filters
---------------

The extraction process marks fragments as low quality ("bogus") if they fail any of the following filters. Bogus fragments are kept in the ``.frag`` file but excluded from downstream analyses.

- **MAPQ ≥ 20**: Minimum mapping quality of every read (``--min-mapq``)
- **Length < 900 bp**: Fragments at or above the bound are excluded (``--max-insert-size``)
- **Valid chromosomes**: Only standard chromosomes (1-22, X, Y, M)
- **No 'N' bases**: No ambiguous nucleotides in the observed fragment end motifs
- **Consistent pairing**: For paired-end data, both reads must map to the same chromosome in opposite orientations with matching template lengths

Duplicate-flagged reads are always skipped; secondary and supplementary alignments only with ``--skip-supplementary``.

General Tips
------------

If you want to process a large set of samples, we recommend against using the ``--bam-dir`` flag. Even though the latter works on *BAM* files in parallel and writes fragment data to disk immediately, it still requires potentially large amounts of memory. We much rather recommend adopting the Nextflow pipeline in our ``tools/`` directory for working with many samples.
