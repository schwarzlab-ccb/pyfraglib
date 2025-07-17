"""
Fragment Processing and BAM File Analysis
=========================================

This module provides the core data structures and functions for processing
cell-free DNA (cfDNA) fragments from BAM files. It includes the fundamental
``Fragment`` class with support for both paired-end and single-end (Nanopore)
sequencing data.

Key Classes
-----------
- :class:`Fragment`: Individual cfDNA fragment with genomic coordinates and
                     properties
- :class:`FragmentList`: Collection of fragments from a single sample
- :class:`FragmentCollection`: Multi-sample fragment collection for batch
                               processing
- :class:`IntervalTable`: Efficient genomic interval management using interval
                          trees (was mostly used in older versions of certain
                          algorithms, but actively used for test cases)

Core Features
-------------
- **BAM File Processing**: Extract fragments from BAM files
- **Mutation Annotation**: Link VCF variants to BAM reads for mutated fragment
                           identification
- **Quality Control**: Filter bogus fragments and apply quality thresholds
- **Batch Processing**: Parallel processing of multiple BAM files
- **Serialization**: Save/load fragment data in compressed .frag format

Fragment Properties
-------------------
Each Fragment contains:

- **Chromosome**: Standardized name of contig from which the fragment
                  originates
- **Genomic Coordinates**: Start and end positions (0-based)
- **Length**: Fragment length in base pairs
- **End Motifs**: 5' and 3' end sequences (k-mers)
- **Mutation Status**: Whether fragment contains variants
- **Quality Flags**: Bogus fragment detection and quality metrics

Example Usage
-------------
Basic fragment extraction and analysis:

.. code-block:: python

    from pyfraglib.fragment import Fragment, FragmentList

    # Extract fragments from BAM file
    fragments = Fragment.from_bam("sample.bam", "variants.vcf")

    # Access fragment properties
    for fragment in fragments:
        print(f"Fragment at {fragment.chrom}:{fragment.start_pos}-"
              f"{fragment.end_pos}")
        print(f"Length: {fragment.length} bp")
        print(f"5' motif: {fragment.end5p}")
        print(f"3' motif: {fragment.end3p}")
        print(f"Mutated: {fragment.is_mutated}")

    # Filter out bogus fragments
    clean_fragments = FragmentList([f for f in fragments if not f.is_bogus])

    # Save to file
    clean_fragments.to_frag_file("sample_clean", "output/")

Batch Processing
----------------
Process multiple BAM files efficiently:

.. code-block:: python

    from pyfraglib.fragment import Fragment

    # Process multiple BAM files in parallel
    bam_files = ["sample1.bam", "sample2.bam", "sample3.bam"]
    vcf_files = ["sample1.vcf", "sample2.vcf", "sample3.vcf"]

    # Batch processing to .frag files.
    Fragment.bams_to_frags(bam_files, vcf_files, "output/")

    # Load all fragments into memory
    fragment_collection = Fragment.from_bams(bam_files, vcf_files)

Importantly, we discourage this method of processing larger sets of samples.
With >20 samples, we recommend using the Nextflow pipeline we are providing.

Interval Operations
-------------------
Efficient genomic interval queries:

.. code-block:: python

    from pyfraglib.fragment import IntervalTable

    # Create interval table
    intervals = IntervalTable()
    intervals.add_bed_file("regions.bed")

    # Find overlapping fragments
    overlapping_fragments = intervals.find_overlapping_fragments(fragments)

    # Check if fragment overlaps any interval
    for fragment in fragments:
        if intervals.overlaps(
            fragment.chrom, fragment.start_pos, fragment.end_pos
        ):
            print(f"Fragment {fragment} overlaps with intervals")

Quality Control
---------------
Fragment quality assessment and filtering:

.. code-block:: python

    from pyfraglib.fragment import Fragment, FragmentList

    # Extract with quality filtering
    fragments = Fragment.from_bam("sample.bam", "variants.vcf")

    # Count quality metrics
    total_fragments = len(fragments)
    bogus_fragments = sum(1 for f in fragments if f.is_bogus)
    mutated_fragments = sum(1 for f in fragments if f.is_mutated)

    print(f"Total fragments: {total_fragments}")
    print(f"Bogus fragments: {bogus_fragments} "
          f"({bogus_fragments/total_fragments:.1%})")
    print(f"Mutated fragments: {mutated_fragments} "
          f"({mutated_fragments/total_fragments:.1%})")

    # Apply additional filtering
    high_quality_fragments = FragmentList([
        f for f in fragments
        if not f.is_bogus and f.length >= 50 and f.length <= 500
    ])

Nanopore Support
----------------
Support for single-end long-read sequencing is implemented generically, but
only tested using Nanopore datasets. Thus, the respective options are named
"nanopore" instead of "long-read, single-ended sequencing":

.. code-block:: python

    from pyfraglib.fragment import Fragment

    # Process Nanopore BAM file
    nanopore_fragments = Fragment.from_bam("nanopore.bam", is_nanopore=True)

Constants
---------
- :const:`INSERT_SIZE_UPPER_BOUND`: Maximum allowed insert size (900 bp)
- :const:`MIN_MAPQ`: Minimum mapping quality threshold (20)
- :const:`DEFAULT_KMER_LEN`: Default k-mer length for end motifs (3)
- :const:`MAX_KMER_LEN`: Maximum k-mer length for end motifs (4)
- :const:`VALID_CHROMOSOME_NAMES`: List of valid chromosome names

Thread Safety
-------------
The Fragment class is immutable and thus thread-safe. Collection classes use
appropriate locking mechanisms for concurrent access. BAM file processing uses
multiprocessing for parallel execution across multiple files - which again is
not recommended for larger datasets (see above).

Performance Considerations
--------------------------
- Use :func:`bams_to_frags` for large-scale processing to avoid memory issues;
  or even better, use the provided Nextflow pipeline
- Fragment collections are memory-intensive
- Interval trees provide O(log n) query performance for genomic overlaps

License
-------
This file is part of ``pyfraglib``, a software suite to calculate fragmentomics
features from cfDNA and perform downstream analyses.

Copyright (C) 2024 Daniel Schütte, daniel.schuette@iccb-cologne.org

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details. You should have received a copy of the GNU General Public
License along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
import gzip
import logging
import os
import pickle
import pysam
import tqdm

from collections import defaultdict
from dataclasses import dataclass
from functools import partial
from intervaltree import IntervalTree, Interval
from multiprocessing import Pool
from pyfraglib.core import fail, PyfragManager, detect_cpus, \
                           homogenize_to_chrom_naming_convention, \
                           get_logger
from typing import Callable, Final, Generator, Optional, Never, cast

INSERT_SIZE_UPPER_BOUND: Final = 900
MIN_MAPQ: Final = 20
DEFAULT_KMER_LEN: Final = 3
MAX_KMER_LEN: Final = 4
VALID_CHROMOSOME_NAMES: Final = [str(x) for x in range(1, 23)] + \
                                ["chr{}".format(x) for x in range(1, 23)] + \
                                ["X", "x", "Y", "y", "chrX", "chrY", "M", "m"]


@dataclass(slots=True)
class Fragment:
    start_pos: int  # first aligned position (0-based)
    end_pos: int  # one past last aligned position
    chrom: str
    is_forward: bool
    length: int
    end5p: str  # 5' -> 3' orientation
    end3p: str  # 5' -> 3' orientation
    is_bogus: bool
    is_mutated: Optional[bool]

    def __init__(
        self, read: pysam.AlignedSegment, mate: Optional[pysam.AlignedSegment],
        mutated_reads: Optional[set[str]]
    ) -> None:
        """
        Initialize a Fragment from BAM alignment data.

        Creates a Fragment object from pysam AlignedSegment data, handling both
        paired-end and single-end sequencing reads. Automatically determines
        fragment coordinates, end motifs, quality flags, and mutation status.

        Args:
            read: Primary aligned read from BAM file. Must be a valid pysam
                AlignedSegment with proper alignment information.
            mate: Mate read for paired-end data, or None for single-end data.
                When None, fragment is processed as single-ended, i.e. Nanopore
            mutated_reads: Set of read names that contain mutations, or None
                if mutation annotation is not available. Used to mark fragments
                as mutated based on VCF analysis.

        Note:
            * Paired-end processing requires both reads to be properly aligned
            * Single-end processing uses the full aligned sequence of the read
            * Fragment coordinates follow BAM 0-based conventions
            * Quality filtering and bogus detection are applied automatically
        """
        if not mate:
            self._init_single_ended(read, mutated_reads)
        else:
            self._init_paired_ended(read, mate, mutated_reads)

    def _init_paired_ended(
        self, read: pysam.AlignedSegment, mate: pysam.AlignedSegment,
        mutated_reads: Optional[set[str]]
    ) -> None:
        # @NOTE(ds): We do all genomic calculations for paired-end reads in one
        # big function. All variables above are filled in after this call.
        self._assign_read_pair_coords(read, mate)

        if mutated_reads is None:
            self.is_mutated = None
        else:
            self.is_mutated = (read.query_name in mutated_reads or
                               mate.query_name in mutated_reads)

        self.is_bogus = self._determine_bogus(read, mate)

    def _init_single_ended(
        self, read: pysam.AlignedSegment,
        mutated_reads: Optional[set[str]]
    ) -> None:
        assert read.reference_end
        assert read.reference_name

        self.start_pos = read.reference_start
        self.end_pos = read.reference_end
        self.chrom = read.reference_name
        self.length = self.end_pos - self.start_pos
        self.is_mutated = (read.query_name in mutated_reads
                           if mutated_reads else None)
        self.is_forward = read.is_forward

        default: str = 'N' * MAX_KMER_LEN
        if read.query_sequence:
            self.end5p = read.query_sequence[0:MAX_KMER_LEN]
            self.end3p = read.query_sequence[-MAX_KMER_LEN:]
        else:
            self.end3p = self.end5p = default

        self.is_bogus = self._determine_bogus(read, None)

    def _determine_bogus(
        self, read: pysam.AlignedSegment, mate: Optional[pysam.AlignedSegment]
    ) -> bool:
        """
        Determine if a fragment should be marked as bogus (low quality).
        Importantly, this method must be called after *all* fragment
        initialization code has run.

        Applies comprehensive quality filters to identify fragments that should
        be excluded from downstream analysis. Filters include mapping quality,
        fragment length, sequence composition, and pairing consistency checks.

        Quality Filters Applied:

        * **Sequence quality**: Fragments with 'N' bases in end motifs
        * **Chromosome validation**: Non-standard chromosome names
        * **Length constraints**: Fragments outside 1-899 bp range
        * **Mapping quality**: MAPQ < 20 for any read
        * **Pairing consistency**: Mismatched orientations or cross-chromosome
                                   pairs
        * **Template length**: Inconsistent insert sizes between mates

        Args:
            read: Primary aligned read for quality assessment.
            mate: Mate read for paired-end validation, or None for single-end.

        Returns:
            bool: True if fragment fails quality filters, False otherwise

        Note:
            * Single-end and paired-end data use different validation criteria
            * Some quality thresholds are defined by module constants
            * Failed fragments are marked but not automatically removed
        """
        is_bogus: bool = False

        if 'N' in self.end5p or 'N' in self.end3p:
            is_bogus = True
        if self.chrom not in VALID_CHROMOSOME_NAMES:
            is_bogus = True

        # @NOTE(ds): We are checking a subset of properties that apply to
        # single-ended data and return early.
        if not mate:
            if self.length >= INSERT_SIZE_UPPER_BOUND or self.length <= 0:
                is_bogus = True
            if read.mapping_quality < MIN_MAPQ:
                is_bogus = True
            return is_bogus

        # @NOTE(ds): We operate on paired-ended sequencing data from this point
        # on.
        read_len: int = abs(read.template_length)
        mate_len: int = abs(mate.template_length)

        # @NOTE(ds): A whole lot of length-related conditions might define a
        # fragment to be bogus.
        if self.length >= INSERT_SIZE_UPPER_BOUND or self.length <= 0:
            is_bogus = True
        elif (self.length != (self.end_pos - self.start_pos)) or \
                (read_len != mate_len):
            is_bogus = True

        if read.mapping_quality < MIN_MAPQ or mate.mapping_quality < MIN_MAPQ:
            is_bogus = True

        # @NOTE(ds): Pairing across contigs and unmatched read orientations
        # are also unexpected and constitute bogus'ness.
        if (read.reference_name != mate.reference_name):
            is_bogus = True
        elif (read.is_read1 and mate.is_read1) or \
                (read.is_read2 and mate.is_read2):
            is_bogus = True

        return is_bogus

    def _assign_read_pair_coords(
        self, read1: pysam.AlignedSegment, read2: pysam.AlignedSegment
    ) -> None:
        # @NOTE(ds): We are only looking at properly aligned, paired reads so
        # we should never encounter the following error.
        if not (read1.reference_end and read2.reference_end and
                read1.reference_name):
            fail("_assign_read_pair_coords: possibly unmapped read")
        s1: int = read1.reference_start
        s2: int = read2.reference_start
        e1: int = read1.reference_end
        e2: int = read2.reference_end
        assert (s1 < e1) and (s2 < e2)

        left_read: pysam.AlignedSegment
        right_read: pysam.AlignedSegment
        if s1 < s2:
            left_read = read1
            right_read = read2
        else:
            # @BUG(ds): Unclear how we should handle cases of overlapping
            # reads (i.e. ``s1 == s2``).
            left_read = read2
            right_read = read1

        self.start_pos = min(s1, s2)
        self.end_pos = max(e1, e2)
        self.chrom = read1.reference_name
        self.length = abs(left_read.template_length)

        # @BUG(ds): Fragile code ahead!
        # @BUG(ds): We did not choose ``read.query_alignment_sequence``. See
        # the ``pysam`` docs for info on differences between the two. Also,
        # sometimes there are Ns in the sequence.
        assert left_read.query_sequence
        assert right_read.query_sequence

        # We could consider different read pair orientations. It seems to be
        # the safest (i.e. most conservative) option to just consider FR
        # read pairs, though.
        default: str = 'N' * MAX_KMER_LEN
        if left_read.is_forward and right_read.is_reverse:
            # According to ``pysam``, the reverse read is already reverse-
            # complemented. Thus, we shouldn't have to do anything, really.
            self.end5p = left_read.query_sequence[0:MAX_KMER_LEN]
            self.end3p = right_read.query_sequence[-MAX_KMER_LEN:]
        else:
            self.end5p = self.end3p = default

        # @BUG(ds): Fragile code ahead!
        if left_read.is_read1:
            self.is_forward = left_read.is_forward
        else:
            self.is_forward = right_read.is_forward

    def dump(self, pickle_file: gzip.GzipFile) -> None:
        pickle.dump(self, pickle_file)

    # @NOTE(ds): In the returned tuple, the first member is 5' and the
    # second member 3'.
    def get_end_motifs(self, kmer_len: int) -> tuple[str, str]:
        """
        Extract k-mer sequences from fragment ends.

        Retrieves nucleotide sequences from the 5' and 3' ends of the fragment
        for motif analysis. These end motifs reflect nuclease cleavage patterns
        and can be used for diversity analysis and tissue-specific signatures.

        Args:
            kmer_len: Length of k-mer sequences to extract. Must be between
                1 and MAX_KMER_LEN (4). Invalid lengths are automatically
                adjusted to DEFAULT_KMER_LEN (3) with a warning.

        Returns:
            tuple[str, str]: A tuple containing (5p_motif, 3p_motif).
                Both motifs are nucleotide sequences of the specified length.

        Note:
            * Motifs are extracted from stored end sequences during fragment
              creation
            * Invalid k-mer lengths trigger automatic correction and logging
            * Used primarily for motif diversity and cleavage pattern analysis

        Example:
            Extract 4-mer end motifs::

                # Get 4-mer sequences from fragment ends
                motif_5p, motif_3p = fragment.get_end_motifs(4)
                print(f"5' motif: {motif_5p}")
                print(f"3' motif: {motif_3p}")
        """
        logger: logging.Logger = get_logger()
        if (kmer_len < 0) or (kmer_len > MAX_KMER_LEN):
            logger.warning(
                f"invalid kmer length {kmer_len}, using {DEFAULT_KMER_LEN}"
            )
            kmer_len = DEFAULT_KMER_LEN
        kmer5p: str = self.end5p[0:kmer_len]
        kmer3p: str = self.end3p[-kmer_len:]
        return (kmer5p, kmer3p)

    @staticmethod
    def from_bam(filepath: str, vcfpath: Optional[str],
                 is_nanopore: bool = False) -> "FragmentList":
        """
        Extract fragments from a BAM file with optional mutation annotation.

        Processes a BAM file to extract cfDNA fragments, applying quality
        filtering and optionally annotating fragments that contain mutations
        from a VCF file. Supports both paired-end and single-end (Nanopore)
        sequencing data.

        Processing includes:
        * **Quality filtering**: MAPQ >= 20, proper fragment lengths
        * **Duplicate handling**: Automatic detection and filtering
        * **Mutation annotation**: Links VCF variants to BAM reads
        * **Progress tracking**: Visual progress bar for large files
        * **Memory optimization**: Efficient read caching and cleanup

        Args:
            filepath: Path to indexed BAM file. Must have corresponding
                .bai index.
            vcfpath: Optional path to VCF file for mutation annotation.
                If provided, fragments containing variants will be marked as
                mutated.
            is_nanopore: Whether to process as single-end Nanopore data.
                Default False assumes paired-end Illumina data.

        Returns:
            FragmentList: List of extracted fragments with quality metrics
                and optional mutation annotations.

        Raises:
            SystemExit: If BAM file lacks required index or processing fails.

        Example:
            Extract fragments with mutation annotation::

                from pyfraglib.fragment import Fragment

                # Basic extraction without mutations
                fragments = Fragment.from_bam("sample.bam", None)

                # With mutation annotation
                fragments = Fragment.from_bam("sample.bam", "variants.vcf")

                print(f"Extracted {len(fragments)} fragments")
                print(f"Bogus fragments: {fragments.count_bogus_fragments()}")

        See Also:
            * :func:`Fragment.from_bams` - Process multiple BAM files in
                                           parallel
            * :func:`Fragment.bams_to_frags` - Batch processing to .frag files
        """
        logger: logging.Logger = get_logger()

        bam_file: pysam.AlignmentFile = pysam.AlignmentFile(filepath)
        if not bam_file.has_index():
            fail("please create an index for BAM file ``{}``".format(filepath))

        mut_reads: Optional[set[str]] = None
        if vcfpath:
            vcf_file: pysam.VariantFile = pysam.VariantFile(vcfpath)
            mut_reads = Fragment.build_mutated_reads_set(bam_file, vcf_file)

        # @NOTE(ds): ``idxstats`` should not fail because a BAI must exist.
        idxstats_output: str = pysam.idxstats(filepath)  # type: ignore
        idxstats_total_reads: int = sum(
            int(line.split('\t')[2]) + int(line.split('\t')[3])
            for line
            in idxstats_output.splitlines()
        )
        # @NOTE(ds): In the following loops, we don't want to do updates for
        # every read, even though ``tqdm`` is supposed to be very fast.
        increment: int = 100_000
        progress_bar: tqdm.tqdm[Never]

        num_total_reads: int = 0
        num_duplicates: int = 0
        fragments: FragmentList = FragmentList()

        if is_nanopore:
            logger.info("processing BAM file ``{}`` as single-ended".format(
                filepath))
            progress_bar = tqdm.tqdm(total=idxstats_total_reads, leave=False)

            for read in bam_file.fetch():
                assert read.query_name

                num_total_reads += 1
                if (num_total_reads % increment == 0):
                    progress_bar.update(increment)

                if read.is_duplicate:
                    num_duplicates += 1
                else:
                    fragments.append(Fragment(read, None, mut_reads))

            progress_bar.close()
            logger.info(
                "total reads: {}, duplicates: {}%".format(
                    num_total_reads,
                    round(100*num_duplicates/num_total_reads, 3)))
        else:
            logger.info("processing BAM file ``{}`` as paired-ended".format(
                filepath))
            progress_bar = tqdm.tqdm(total=idxstats_total_reads, leave=False)

            num_unmapped: int = 0
            read_cache: dict[str, pysam.AlignedSegment] = {}

            for read in bam_file.fetch():
                assert read.query_name

                invalid_read: bool = False
                num_total_reads += 1

                if (num_total_reads % increment == 0):
                    progress_bar.update(increment)

                # @NOTE(ds): We cannot use ``is_proper_pair`` here because some
                # aligners do weird stuff and set flag 0x2 incorrectly. We
                # just find our pairs ourselves.
                if read.is_unmapped:
                    num_unmapped += 1
                    invalid_read = True
                if read.is_duplicate:
                    num_duplicates += 1
                    invalid_read = True

                if not invalid_read:
                    if read.query_name in read_cache:
                        mate: pysam.AlignedSegment = \
                            read_cache[read.query_name]
                        fragments.append(Fragment(read, mate, mut_reads))
                        read_cache.pop(read.query_name)  # _major_ mem saver
                    else:
                        read_cache[read.query_name] = read

            progress_bar.close()
            logger.info(
                "total reads: {}, unmapped: {:.3}%, duplicates: {:.3}%".format(
                    num_total_reads, 100*num_unmapped/num_total_reads,
                    100*num_duplicates/num_total_reads))
            logger.info(
                "{:.3}% reads left in cache (probably unpaired)".format(
                    100*len(read_cache)/num_total_reads))

        # @NOTE(ds): Careful cleanup. We must not leak memory.
        bam_file.close()
        if vcfpath:
            vcf_file.close()

        return fragments

    @staticmethod
    def build_mutated_reads_set(
        bam_file: pysam.AlignmentFile, vcf_file: pysam.VariantFile
    ) -> set[str]:
        """
        Build a set of read names that contain mutations from VCF analysis.

        Links VCF variant records to BAM reads by checking which reads contain
        alternative alleles at variant positions. This enables identification
        of cfDNA fragments carrying specific mutations for downstream analysis.

        The function:
        * Processes all SNVs in the VCF file
        * Finds overlapping reads at each variant position
        * Checks read bases against reference and alternative alleles
        * Returns read names that contain mutations

        Args:
            bam_file: Opened BAM file containing aligned reads. Must be indexed
                and coordinate-sorted for efficient variant lookups.
            vcf_file: Opened VCF file containing variant calls. Only SNVs
                (single nucleotide variants) are processed.

        Returns:
            set[str]: Set of read names that contain alternative alleles at
                variant positions. Used to mark fragments as mutated.

        Note:
            * Only processes single nucleotide variants (SNVs)
            * Handles chromosome naming inconsistencies automatically
            * Tracks duplex sequencing support for quality assessment
              (very specific to our workflow, thus not of interest to most
              other researchers)
            * Skips reads with ambiguous bases (neither ref nor alt)

        Example:
            Build mutated reads set for fragment annotation::

                import pysam

                # Open files
                bam_file = pysam.AlignmentFile("sample.bam")
                vcf_file = pysam.VariantFile("variants.vcf")

                # Build mutated reads set
                mutated_reads = Fragment.build_mutated_reads_set(
                    bam_file, vcf_file
                )
                print(f"Found {len(mutated_reads)} mutated reads")

                # Close files
                bam_file.close()
                vcf_file.close()
        """
        logger: logging.Logger = get_logger()
        mutated_reads: set[str] = set()

        num_unknown_variants: int = 0
        num_singles: int = 0
        num_mutated_reads: int = 0
        variant: pysam.VariantRecord
        fname: bytes = vcf_file.filename
        vcf_filename: str = fname.decode()
        for variant in vcf_file:
            if variant.rlen != 1:
                logger.warn("skipping variant record of non-SNV (len={}) "
                            "in ``{}``".format(variant.rlen, vcf_filename))
                continue

            # @NOTE(ds): The documentation of ``pysam`` isn't really clear
            # about 0- and 1-based indexing. We try to be consistent and always
            # use the methods that are 0-based.
            read: pysam.AlignedSegment
            contig: str = homogenize_to_chrom_naming_convention(
                variant.contig, bam_file.header.to_dict()  # type: ignore
            )

            try:
                genomic_position: pysam.IteratorRow = bam_file.fetch(
                    contig=contig, start=variant.start, stop=variant.stop
                )
            except ValueError as e:
                fail(f"failed to fetch {variant} from bam file: {e}")

            for read in genomic_position:
                read_positions: list[int] = \
                    read.get_reference_positions(full_length=True)

                if variant.start in read_positions:
                    read_index: int = read_positions.index(variant.start)

                    assert read.query_sequence
                    assert variant.alts

                    read_base: str = read.query_sequence[read_index]
                    if read_base in variant.alts:
                        if not is_duplex(read):
                            num_singles += 1
                        num_mutated_reads += 1
                        if read.query_name:
                            mutated_reads.add(read.query_name)
                    elif read_base != variant.ref:
                        # NOTE(ds): Read has a base that's neither a ref nor
                        # an alt allele. We _do not_ include this read in the
                        # mutated set because the mutation is not well-defined.
                        # It was probably filtered out during variant calling.
                        num_unknown_variants += 1
                else:
                    # The variant position is probably deleted from ``read``.
                    pass

        logger.info("skipped at least {} variants because read base did not "
                    "match ref nor alt allele (``{}``)".format(
                        num_unknown_variants, vcf_filename))
        logger.info("found {} mutated reads without duplex support (out of "
                    "{} mutated reads total)".format(num_singles,
                                                     num_mutated_reads))

        return mutated_reads

    @staticmethod
    def from_bams(
        filepaths: list[str], vcfpaths: Optional[list[str]],
        is_nanopore: bool = False
    ) -> "FragmentCollection":
        """
        Process multiple BAM files in parallel and collect fragments in memory.

        Warning:
            This function has significant memory overhead. Consider using
            :func:`Fragment.bams_to_frags` for large-scale processing or the
            provided Nextflow pipeline.
        """
        with PyfragManager() as mngr:
            shared_collection = mngr.FragmentCollection()  # type: ignore
            input_data: list[tuple[str, str | None]] = []

            for idx, filepath in enumerate(filepaths):
                vcfpath: Optional[str] = \
                    None if not vcfpaths else vcfpaths[idx]
                input_data.append((filepath, vcfpath))

            with Pool(processes=detect_cpus()) as pool:
                partial_task: Callable[[str, str], None] = partial(
                    task0, frags_per_bam=shared_collection,  # type: ignore
                    is_nanopore=is_nanopore
                )
                pool.starmap(partial_task, input_data)

            return cast(
                FragmentCollection,
                shared_collection._getvalue()  # type: ignore
            )

    @staticmethod
    def bams_to_frags(
        filepaths: list[str], vcfpaths: Optional[list[str]], out_dir: str,
        is_nanopore: bool = False
    ) -> None:
        """
        Process multiple BAM files in parallel and save fragments to .frag
        files.

        Batch processing that extracts fragments from multiple BAM files in
        parallel and writes results directly to disk as compressed .frag files.
        This approach avoids the memory overhead of :func:`Fragment.from_bams`
        by not storing all fragments in memory simultaneously. It can still be
        very expensive.

        Note:
            We use the same multiprocessing idioms as in ``from_bams``
        """
        input_data: list[tuple[str, str | None]] = []

        for idx, filepath in enumerate(filepaths):
            vcfpath: Optional[str] = \
                None if not vcfpaths else vcfpaths[idx]
            input_data.append((filepath, vcfpath))
        with Pool(processes=detect_cpus()) as pool:
            partial_task: Callable[[str, str], None] = partial(
                task1, out_dir=out_dir, is_nanopore=is_nanopore
            )
            pool.starmap(partial_task, input_data)

    @classmethod
    def create_simulated(
        cls, start_pos: int, end_pos: int, chrom: str, length: int,
        end5p: str, end3p: str, is_forward: bool = True,
        is_mutated: bool | None = None
    ) -> "Fragment":
        """
        Create a simulated Fragment without requiring pysam objects.

        Creates a Fragment object with known properties for simulation
        purposes, bypassing the normal BAM file processing pipeline. This
        method is primarily used by the simulation module to generate synthetic
        cfDNA fragments with realistic characteristics.

        Args:
            start_pos: Fragment start position (0-based genomic coordinate).
            end_pos: Fragment end position (0-based, exclusive).
            chrom: Chromosome name (e.g., "chr1", "1").
            length: Fragment length in base pairs.
            end5p: 5' end motif sequence (nucleotide string).
            end3p: 3' end motif sequence (nucleotide string).
            is_forward: Strand orientation. Default True for forward strand.
            is_mutated: Mutation status, or None if unknown.

        Returns:
            Fragment: A properly initialized Fragment object with simulated
                properties.

        Note:
            * Automatically applies quality filtering (bogus detection)
            * Bypasses BAM file dependency for simulation workflows
            * Used primarily by :mod:`pyfraglib.simulator` modules
        """
        fragment: Fragment = cls.__new__(cls)
        fragment.start_pos = start_pos
        fragment.end_pos = end_pos
        fragment.chrom = chrom
        fragment.length = length
        fragment.end5p = end5p
        fragment.end3p = end3p
        fragment.is_forward = is_forward
        fragment.is_mutated = is_mutated

        fragment.is_bogus = (
            "N" in end5p or "N" in end3p or
            chrom not in VALID_CHROMOSOME_NAMES or
            length <= 0 or length >= INSERT_SIZE_UPPER_BOUND
        )

        return fragment


# @NOTE(ds): Constraints of multiprocessing and pickle force us to define these
# functions outside of ``from_bams`` and ``bams_to_frags``. Tasks 0 and 1 only
# differ in how they treat the resulting ``FragmentList``s. Whereas t0 writes
# the into a shared buffer, t1 writes them to disk immediately. That's freeing
# the used memory and reducing our mem footprint.
def task0(filepath: str, vcfpath: str | None,
          frags_per_bam: "FragmentCollection", is_nanopore: bool) -> None:
    frags: FragmentList = Fragment.from_bam(filepath, vcfpath, is_nanopore)
    filename: str = os.path.basename(filepath)
    name, _ = os.path.splitext(filename)
    frags_per_bam.append(name, frags)


def task1(filepath: str, vcfpath: str | None,
          out_dir: str, is_nanopore: bool) -> None:
    logger: logging.Logger = get_logger()

    frags: FragmentList = Fragment.from_bam(filepath, vcfpath, is_nanopore)
    filename: str = os.path.basename(filepath)
    name, _ = os.path.splitext(filename)

    # @NOTE(ds): We are pretty inconsistent when it comes to logging. Sometimes
    # we log in the worker functions (like we do in ``from_bam``), but
    # sometimes we log in the dispatchers, too (like below).
    logger.info("saving ``{}.frag`` to ``{}``".format(name, out_dir))
    frags.to_frag_file(name, out_dir)


@dataclass(slots=True)
class FragmentList():
    """
    Collection of Fragment objects from a single sample with analysis methods.

    A container class for managing collections of Fragment objects extracted
    from a single BAM file. Provides methods for fragment manipulation,
    analysis, and serialization to .frag files for efficient storage and reuse.

    Key Features:
    * **Fragment management**: Add, iterate, and count fragments
    * **Quality metrics**: Count bogus and mutated fragments
    * **Motif analysis**: Extract and count end motif patterns
    * **Serialization**: Save/load to compressed .frag format
    * **Interval operations**: Convert to interval trees for overlap queries

    Example:
        Basic fragment collection usage::

            from pyfraglib.fragment import Fragment, FragmentList

            # Create empty collection
            fragments = FragmentList()

            # Add fragments (typically from Fragment.from_bam)
            fragments.append(fragment1)
            fragments.append(fragment2)

            # Get statistics
            total = fragments.length()
            bogus = fragments.count_bogus_fragments()
            mutated = fragments.count_mutated_fragments()

            print(f"Total: {total}, Bogus: {bogus}, Mutated: {mutated}")

            # Save to file
            fragments.to_frag_file("sample", "output/")

    See Also:
        * :class:`Fragment` - Individual fragment objects
        * :class:`FragmentCollection` - Multi-sample fragment collections
        * :class:`FragFile` - Reading .frag files back into memory
    """
    __fragments: list[Fragment]

    def __init__(self) -> None:
        self.__fragments: list[Fragment] = []

    def append(self, fragment: Fragment) -> None:
        self.__fragments.append(fragment)

    def length(self) -> int:
        return len(self.__fragments)

    def __iter__(self) -> Generator[Fragment, None, None]:
        for fragment in self.__fragments:
            yield fragment

    def count_bogus_fragments(self) -> int:
        counter: int = 0
        for frag in self.__fragments:
            if frag.is_bogus:
                counter += 1
        return counter

    def count_mutated_fragments(self) -> int:
        counter: int = 0
        for frag in self.__fragments:
            if frag.is_mutated:
                counter += 1
        return counter

    def to_frag_file(self, name: str, out_dir: str) -> None:
        outfile_path: str = os.path.join(out_dir, name) + ".frag"
        outfile: gzip.GzipFile
        with gzip.open(outfile_path, "wb") as outfile:
            for fragment in self.__fragments:
                fragment.dump(outfile)

    # @NOTE(ds): Returns the 5' and 3' motifs as well as the number of
    # analyzed fragments.
    def count_endmotifs(
        self, kmer_len: int
    ) -> tuple[defaultdict[str, int], defaultdict[str, int], int]:
        motifs_5p: defaultdict[str, int] = defaultdict(int)
        motifs_3p: defaultdict[str, int] = defaultdict(int)

        num_frags: int = 0
        fragment: Fragment

        for fragment in self.__fragments:
            if fragment.is_bogus:
                continue

            num_frags += 1
            end5p, end3p = fragment.get_end_motifs(kmer_len)
            motifs_5p[end5p] += 1
            _ = motifs_3p[end5p]
            motifs_3p[end3p] += 1
            _ = motifs_5p[end3p]

        return motifs_5p, motifs_3p, num_frags

    def to_interval_table(self) -> "IntervalTable":
        table: IntervalTable = IntervalTable()

        for frag in self.__fragments:
            if frag.is_bogus:
                continue
            table.insert(frag)

        return table


def is_duplex(read: pysam.AlignedSegment) -> bool:
    has_xi: bool = True
    has_xj: bool = True

    try:
        read.get_tag("XI")
    except KeyError:
        has_xi = False

    try:
        read.get_tag("XJ")
    except KeyError:
        has_xj = False

    return has_xi and has_xj


class IntervalTable():
    """
    Efficient genomic interval management using chromosome-indexed interval
    trees.

    Provides fast genomic interval operations by maintaining separate interval
    trees for each chromosome. Enables efficient overlap queries and
    intersection detection.

    Key Features:
    * **Chromosome indexing**: Separate interval trees per chromosome
    * **Efficient queries**: O(log n) overlap detection performance
    * **Fragment integration**: Direct fragment insertion and querying

    Example:
        Basic interval operations::

            from pyfraglib.fragment import IntervalTable, Fragment

            # Create interval table
            intervals = IntervalTable()

            # Insert fragments
            fragments = Fragment.from_bam("sample.bam")
            for fragment in fragments:
                if not fragment.is_bogus:
                    intervals.insert(fragment)

            # Query overlaps in a genomic window
            overlaps = intervals.get_overlaps("chr1", 1000000, 1001000)
            print(f"Found {len(overlaps)} overlapping intervals")

    Note:
        * Interval trees are created on-demand for each chromosome
        * Uses 0-based coordinates following BAM conventions

    See Also:
        * :meth:`FragmentList.to_interval_table` - Convert fragments to
                                                   intervals
    """
    def __init__(self) -> None:
        self.__table: dict[str, IntervalTree] = {}  # type: ignore

    def insert(self, fragment: Fragment) -> None:
        if fragment.chrom not in self.__table:  # type: ignore
            self.__table[fragment.chrom] = IntervalTree()  # type: ignore

        self.__table[fragment.chrom].addi(  # type: ignore
            fragment.start_pos, fragment.end_pos
        )

    def get_overlaps(  # type: ignore
        self, chrom: str, win_start: int, win_end: int
    ) -> list[Interval]:
        itree: IntervalTree = self.__table[chrom]  # type: ignore
        intervals: list[Interval] = itree[win_start:win_end]  # type: ignore
        return intervals  # type: ignore


class FragmentCollection():
    """
    Multi-sample container for FragmentList collections with batch operations.

    Dictionary-like container that manages FragmentList objects from multiple
    samples, enabling cross-sample analysis and batch processing operations.
    Typically used for studies involving multiple BAM files or comparative
    analysis workflows. Obviously very memory-intensive as soon as the number
    of samples grows larger.

    Key Features:
    * **Multi-sample management**: Store fragments from multiple samples
    * **Dictionary interface**: Access samples by name with iteration support
    * **Batch operations**: Apply operations across all samples

    Example:
        Multi-sample fragment analysis::

            from pyfraglib.fragment import Fragment, FragmentCollection

            # Load multiple samples
            bam_files = ["sample1.bam", "sample2.bam", "sample3.bam"]
            collection = Fragment.from_bams(bam_files, None)

            # Iterate through samples
            for sample_name, fragments in collection:
                print(f"Sample {sample_name}: {len(fragments)} fragments")
                print(f"Bogus: {fragments.count_bogus_fragments()}")
                print(f"Mutated: {fragments.count_mutated_fragments()}")

            # Batch save to .frag files
            collection.to_frag_files("output/")

            # Access sample names
            sample_names = collection.record_names()
            print(f"Loaded {len(sample_names)} samples: {sample_names}")

    Warning:
        Large collections can consume significant memory. Consider using
        :func:`Fragment.bams_to_frags` for memory-efficient batch processing
        or the Nextflow pipeline.

    See Also:
        * :func:`Fragment.from_bams` - Create collections from multiple BAMs
        * :func:`Fragment.bams_to_frags` - Memory-efficient batch processing
        * :class:`FragmentList` - Individual sample fragment collections
    """
    def __init__(self) -> None:
        self.__fragments_by_name: dict[str, FragmentList] = dict()

    def append(self, name: str, fragments: FragmentList) -> None:
        self.__fragments_by_name[name] = fragments

    def __iter__(self) -> Generator[tuple[str, FragmentList], None, None]:
        for name, fragments in self.__fragments_by_name.items():
            yield (name, fragments)

    def record_names(self) -> list[str]:
        return list(self.__fragments_by_name.keys())

    def to_frag_files(self, out_dir: str) -> None:
        for name, fragments in self:
            fragments.to_frag_file(name, out_dir)
