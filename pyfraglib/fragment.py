# This file is part of `pyfraglib`, a software suite to calculate fragmentomics
# features from cfDNA and perform downstream analyses.
#
# Copyright (C) 2024 Daniel Schütte, daniel.schuette@iccb-cologne.org
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details. You should have received a copy of the GNU General Public
# License along with this program. If not, see <https://www.gnu.org/licenses/>.
import gzip
import logging
import os
import pickle
import pysam
import tqdm

from collections import defaultdict
from dataclasses import dataclass
from functools import partial
from intervaltree import IntervalTree, Interval  # type: ignore
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
    # @NOTE(ds): `left' and `right' are defined according to genomic
    # orientation. We assume that the forward (+) strand (5' -> 3') has
    # increasing coordinates (as it is the case with BAM files). Also,
    # positions are 0-based.
    #
    # @NOTE(ds): We cannot store `pysam.AlignedSegment's in this class
    # because otherwise serialization would not work.
    start_pos: int  # first aligned position (0-based)
    end_pos: int  # one past last aligned position
    chrom: str
    is_forward: bool
    length: int
    end5p: str  # 5' -> 3' orientation
    end3p: str  # 5' -> 3' orientation
    is_bogus: bool
    is_mutated: Optional[bool]

    # @NOTE(ds): With `mate' set to `None', we are assuming single-ended data.
    def __init__(
        self, read: pysam.AlignedSegment, mate: Optional[pysam.AlignedSegment],
        mutated_reads: Optional[set[pysam.AlignedSegment]]
    ) -> None:
        if not mate:
            self.init_single_ended(read, mutated_reads)
        else:
            self.init_paired_ended(read, mate, mutated_reads)

    def init_paired_ended(
        self, read: pysam.AlignedSegment, mate: pysam.AlignedSegment,
        mutated_reads: Optional[set[pysam.AlignedSegment]]
    ) -> None:
        # @NOTE(ds): We do all genomic calculations for paired-end reads in one
        # big function. All variables above are filled in after this call.
        self.assign_read_pair_coords(read, mate)

        if not mutated_reads:
            self.is_mutated = None
        else:
            self.is_mutated = read in mutated_reads or mate in mutated_reads

        self.is_bogus = self.determine_bogus(read, mate)  # must be called last

    def init_single_ended(
        self, read: pysam.AlignedSegment,
        mutated_reads: Optional[set[pysam.AlignedSegment]]
    ) -> None:
        assert read.reference_end
        assert read.reference_name

        self.start_pos = read.reference_start
        self.end_pos = read.reference_end
        self.chrom = read.reference_name
        self.length = self.end_pos - self.start_pos
        self.is_mutated = read in mutated_reads if mutated_reads else None
        self.is_forward = read.is_forward

        default: str = 'N' * MAX_KMER_LEN
        if read.query_sequence:
            self.end5p = read.query_sequence[0:MAX_KMER_LEN]
            self.end3p = read.query_sequence[-MAX_KMER_LEN:]
        else:
            self.end3p = self.end5p = default

        self.is_bogus = self.determine_bogus(read, None)  # must be called last

    def determine_bogus(
        self, read: pysam.AlignedSegment, mate: Optional[pysam.AlignedSegment]
    ) -> bool:
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

    def assign_read_pair_coords(
        self, read1: pysam.AlignedSegment, read2: pysam.AlignedSegment
    ) -> None:
        # @NOTE(ds): We are only looking at properly aligned, paired reads so
        # we should never encounter the following error.
        if not (read1.reference_end and read2.reference_end and
                read1.reference_name):
            fail("assign_read_pair_coords: possibly unmapped read")
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
            # reads (i.e. `s1 == s2').
            left_read = read2
            right_read = read1

        self.start_pos = min(s1, s2)
        self.end_pos = max(e1, e2)
        self.chrom = read1.reference_name
        self.length = abs(left_read.template_length)

        # @BUG(ds): Fragile code ahead!
        # @BUG(ds): We did not choose `read.query_alignment_sequence'. See the
        # `pysam' docs for info on differences between the two. Also, sometimes
        # there are Ns in the sequence.
        assert left_read.query_sequence
        assert right_read.query_sequence

        # We could consider different read pair orientations. It seems to be
        # the savest (i.e. most conservative) option to just consider FR
        # read pairs, though.
        default: str = 'N' * MAX_KMER_LEN
        if left_read.is_forward and right_read.is_reverse:
            # According to `pysam', the reverse read is already reverse-
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
        logger: logging.Logger = get_logger()
        if (kmer_len < 0) or (kmer_len > MAX_KMER_LEN):
            logger.warning(
                f"invalid kmer length {kmer_len}, using {DEFAULT_KMER_LEN}"
            )
            kmer_len = DEFAULT_KMER_LEN
        kmer5p: str = self.end5p[0:kmer_len]
        kmer3p: str = self.end3p[-kmer_len:]
        return (kmer5p, kmer3p)

    # `from_bam' queries a single BAM file for the fragments that it contains.
    # It _does not_ ensure that `filepath' is valid - the caller must do so.
    # The same is true for the optional VCF file that can be provided to enable
    # annotation of mutated fragments.
    @staticmethod
    def from_bam(filepath: str, vcfpath: Optional[str],
                 is_nanopore: bool = False) -> "FragmentList":
        logger: logging.Logger = get_logger()

        bam_file: pysam.AlignmentFile = pysam.AlignmentFile(filepath)
        if not bam_file.has_index():
            fail("please create an index for BAM file `{}'".format(filepath))

        mut_reads: Optional[set[pysam.AlignedSegment]] = None
        if vcfpath:
            vcf_file: pysam.VariantFile = pysam.VariantFile(vcfpath)
            mut_reads = Fragment.build_mutated_reads_set(bam_file, vcf_file)

        # @NOTE(ds): `idxstats' should not fail because a BAI must exist.
        idxstats_output: str = pysam.idxstats(filepath)  # type: ignore
        idxstats_total_reads: int = sum(
            int(line.split('\t')[2]) + int(line.split('\t')[3])
            for line
            in idxstats_output.splitlines()
        )
        # @NOTE(ds): In the following loops, we don't want to do updates every
        # read, even though `tqdm' is supposed to be very fast.
        increment: int = 100_000
        progress_bar: tqdm.tqdm[Never]

        num_total_reads: int = 0
        num_duplicates: int = 0
        fragments: FragmentList = FragmentList()

        if is_nanopore:
            logger.info("processing BAM file `{}' as single-ended".format(
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
            logger.info("processing BAM file `{}' as paired-ended".format(
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

                # @NOTE(ds): We cannot use `is_proper_pair' here because some
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
    ) -> set[pysam.AlignedSegment]:
        logger: logging.Logger = get_logger()
        mutated_reads: set[pysam.AlignedSegment] = set()

        num_unknown_variants: int = 0
        num_singles: int = 0
        num_mutated_reads: int = 0
        variant: pysam.VariantRecord
        fname: bytes = vcf_file.filename
        vcf_filename: str = fname.decode()
        for variant in vcf_file:
            if variant.rlen != 1:
                logger.warn("skipping variant record of non-SNV (len={}) "
                            "in `{}'".format(variant.rlen, vcf_filename))
                continue

            # @NOTE(ds): The documentation of `pysam' isn't really clear about
            # 0- and 1-based indexing. We try to be consistent and always use
            # the methods that are 0-based.
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
                        mutated_reads.add(read)
                    elif read_base != variant.ref:
                        # NOTE(ds): Read has a base that's neither a ref nor
                        # an alt allele. We _do not_ include this read in the
                        # mutated set because the mutation is not well-defined.
                        # It was probably filtered out during variant calling.
                        num_unknown_variants += 1
                else:
                    # The variant position is probably deleted from `read'.
                    pass

        logger.info("skipped at least {} variants because read base did not "
                    "match ref nor alt allele (`{}')".format(
                        num_unknown_variants, vcf_filename))
        logger.info("found {} mutated reads without duplex support (out of "
                    "{} mutated reads total)".format(num_singles,
                                                     num_mutated_reads))

        return mutated_reads

    # `from_bams' builds ontop of `from_bam' by retrieving fragment information
    # from multiple BAM files and storing them in a dictionary. Again, we do
    # not ensure that `filepaths' contains just BAM file paths. Dict keys are
    # just the file names (without the .bam extension and without any directory
    # prefixes). To be a little more sophisticated, we do not return a naked
    # dict but a thin wrapper object.
    # This function comes with a lot of memory overhead! It is not advisable to
    # use it for larger collections of BAM files. The overhead is because
    # individual BAM files are processed and the results are accumulated in
    # memory. Unfortunately, the in-memory representation of `Fragment's is
    # quite large.
    @staticmethod
    def from_bams(
        filepaths: list[str], vcfpaths: Optional[list[str]],
        is_nanopore: bool = False
    ) -> "FragmentCollection":
        # @NOTE(ds): Unfortunately, we have to do a lot of dynamic typing.
        # Thus, mypy complains on multiple occasions.
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

    # @NOTE(ds): As opposed to `from_bams', `bams_to_frags' works on multiple
    # BAM files in parallel, writing them out to FRAG files _without_
    # collecting all data in memory. That's much more efficient. Also, the
    # caller can specify that the input BAM files contain unpair Nanopore
    # reads.
    @staticmethod
    def bams_to_frags(
        filepaths: list[str], vcfpaths: Optional[list[str]], out_dir: str,
        is_nanopore: bool = False
    ) -> None:
        # @NOTE(ds): We use the same multiprocessing idioms as in `from_bams'.
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


# @NOTE(ds): Constraints of multiprocessing and pickle force us to define these
# functions outside of `from_bams' and `bams_to_frags'. Tasks 0 and 1 only
# differ in how they treat the resulting `FragmentList's. Whereas t0 writes
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
    # we log in the worker functions (like we do in `from_bam'), but sometimes
    # we log in the dispatchers, too (like below).
    logger.info("saving `{}.frag' to `{}'".format(name, out_dir))
    frags.to_frag_file(name, out_dir)


@dataclass(slots=True)
class FragmentList():
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


# @NOTE(ds): An `IntervalTable' maps chromosomes to a unique `IntervalTree'.
class IntervalTable():
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
