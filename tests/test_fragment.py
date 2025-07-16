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
import unittest
import tempfile
import os

from unittest.mock import Mock
from typing import Final
from pyfraglib.fragment import Fragment, FragmentList, FragmentCollection, \
                               IntervalTable, is_duplex, \
                               VALID_CHROMOSOME_NAMES, \
                               INSERT_SIZE_UPPER_BOUND, MIN_MAPQ, \
                               DEFAULT_KMER_LEN, MAX_KMER_LEN
from tests.test_fixtures import MockAlignedSegment, create_mock_fragment, \
                                create_test_fragment_list, \
                                create_mock_variant_record, \
                                create_temp_frag_file, cleanup_temp_file

# @NOTE(ds): These files won't be available for anyone but the main developer.
# If the test runner does not find them, the respective integration test is
# ignored.
MUTATED_READS_ONLY_TEST_BAM: Final[str] = \
    "/home/daniel/lab/code/pyfraglib/data/DED005_mutated_only.bam"
MUTATED_READS_ONLY_TEST_BAM_BAI: Final[str] = \
    "/home/daniel/lab/code/pyfraglib/data/DED005_mutated_only.bam.bai"
MUTATED_READS_ONLY_TEST_VCF: Final[str] = \
    "/home/daniel/lab/code/pyfraglib/data/full/DED005_BL_full.vcf"
UNMUTATED_READS_ONLY_TEST_BAM: Final[str] = \
    "/home/daniel/lab/code/pyfraglib/data/DED005_unmutated_only_0.01.bam"
UNMUTATED_READS_ONLY_TEST_BAM_BAI: Final[str] = \
    "/home/daniel/lab/code/pyfraglib/data/DED005_unmutated_only_0.01.bam.bai"
UNMUTATED_READS_ONLY_TEST_VCF: Final[str] = \
    MUTATED_READS_ONLY_TEST_VCF
MIXED_MUTATED_READS_TEST_BAM: Final[str] = \
    "/home/daniel/lab/code/pyfraglib/data/DED005_BL_0.01.bam"
MIXED_MUTATED_READS_TEST_BAM_BAI: Final[str] = \
    "/home/daniel/lab/code/pyfraglib/data/DED005_BL_0.01.bam.bai"
MIXED_MUTATED_READS_TEST_VCF: Final[str] = \
    MUTATED_READS_ONLY_TEST_VCF


class TestFragment(unittest.TestCase):
    """Test Fragment class functionality."""

    def test_fragment_init_single_ended(self) -> None:
        """Test Fragment initialization with single-ended read."""
        read = MockAlignedSegment(
            reference_start=100,
            reference_end=250,
            reference_name="1",
            query_sequence="ATCGATCGATCGATCG",
            mapping_quality=30
        )

        fragment: Fragment = Fragment(read, None, None)  # type: ignore
        self.assertEqual(fragment.start_pos, 100)
        self.assertEqual(fragment.end_pos, 250)
        self.assertEqual(fragment.chrom, "1")
        self.assertEqual(fragment.length, 150)
        self.assertEqual(fragment.end5p, "ATCG")
        self.assertEqual(fragment.end3p, "ATCG")
        self.assertTrue(fragment.is_forward)
        self.assertIsNone(fragment.is_mutated)
        self.assertFalse(fragment.is_bogus)

    def test_fragment_init_paired_ended(self) -> None:
        """Test Fragment initialization with paired-ended reads."""
        read1 = MockAlignedSegment(
            reference_start=100,
            reference_end=200,
            reference_name="1",
            query_sequence="ATCGATCGATCGATCG",
            is_forward=True,
            is_read1=True,
            is_read2=False,
            template_length=150
        )

        read2 = MockAlignedSegment(
            reference_start=175,
            reference_end=250,
            reference_name="1",
            query_sequence="CGATCGATCGATCGAT",
            is_forward=False,
            is_reverse=True,
            is_read1=False,
            is_read2=True,
            template_length=150
        )

        fragment: Fragment = Fragment(read1, read2, None)  # type: ignore
        self.assertEqual(fragment.start_pos, 100)
        self.assertEqual(fragment.end_pos, 250)
        self.assertEqual(fragment.chrom, "1")
        self.assertEqual(fragment.length, 150)
        self.assertEqual(fragment.end5p, "ATCG")
        self.assertEqual(fragment.end3p, "CGAT")
        self.assertIsNone(fragment.is_mutated)
        self.assertFalse(fragment.is_bogus)

    def test_fragment_bogus_detection(self) -> None:
        """Test bogus fragment detection."""
        read_with_n = MockAlignedSegment(
            query_sequence="ATCNGATCGATCGATC"
        )
        fragment: Fragment = Fragment(read_with_n, None, None)  # type: ignore
        self.assertTrue(fragment.is_bogus)

        read_invalid_chr = MockAlignedSegment(
            reference_name="invalid_chr"
        )
        fragment = Fragment(read_invalid_chr, None, None)  # type: ignore
        self.assertTrue(fragment.is_bogus)

        read_low_mapq = MockAlignedSegment(
            mapping_quality=10
        )
        fragment = Fragment(read_low_mapq, None, None)  # type: ignore
        self.assertTrue(fragment.is_bogus)

        read_long = MockAlignedSegment(
            reference_start=100,
            reference_end=1100
        )
        fragment = Fragment(read_long, None, None)  # type: ignore
        self.assertTrue(fragment.is_bogus)

    def test_fragment_mutated_annotation(self) -> None:
        """Test fragment mutation annotation."""
        read = MockAlignedSegment()
        mutated_reads: set[str] = {read.query_name}

        fragment: Fragment = Fragment(
            read, None, mutated_reads  # type: ignore
        )
        self.assertTrue(fragment.is_mutated)

        other_read = MockAlignedSegment(query_name="other_read")
        fragment = Fragment(other_read, None, mutated_reads)  # type: ignore
        self.assertFalse(fragment.is_mutated)

    def test_fragment_get_end_motifs(self) -> None:
        """Test end motif extraction."""
        fragment = create_mock_fragment(
            end5p="ATCGATCG",
            end3p="CGATCGAT"
        )

        motif5p, motif3p = fragment.get_end_motifs(3)
        self.assertEqual(motif5p, "ATC")
        self.assertEqual(motif3p, "GAT")

        motif5p, motif3p = fragment.get_end_motifs(4)
        self.assertEqual(motif5p, "ATCG")
        self.assertEqual(motif3p, "CGAT")

        with self.assertLogs("pyfraglib", level="WARNING"):
            motif5p, motif3p = fragment.get_end_motifs(10)
            expected_len = DEFAULT_KMER_LEN
            self.assertEqual(len(motif5p), expected_len)
            self.assertEqual(len(motif3p), expected_len)


class TestFragmentList(unittest.TestCase):
    """Test FragmentList class functionality."""

    def test_fragment_list_basic_operations(self) -> None:
        """Test basic FragmentList operations."""
        fragment_list = FragmentList()
        self.assertEqual(fragment_list.length(), 0)

        fragment = create_mock_fragment()
        fragment_list.append(fragment)
        self.assertEqual(fragment_list.length(), 1)

        fragments = list(fragment_list)
        self.assertEqual(len(fragments), 1)
        self.assertEqual(fragments[0], fragment)

    def test_fragment_list_count_bogus(self) -> None:
        """Test counting bogus fragments."""
        fragment_list = FragmentList()

        normal_fragment = create_mock_fragment(is_bogus=False)
        fragment_list.append(normal_fragment)
        bogus_fragment = create_mock_fragment(is_bogus=True)
        fragment_list.append(bogus_fragment)

        self.assertEqual(fragment_list.count_bogus_fragments(), 1)

    def test_fragment_list_count_mutated(self) -> None:
        """Test counting mutated fragments."""
        fragment_list = FragmentList()

        normal_fragment = create_mock_fragment(is_mutated=False)
        fragment_list.append(normal_fragment)
        mutated_fragment = create_mock_fragment(is_mutated=True)
        fragment_list.append(mutated_fragment)
        unknown_fragment = create_mock_fragment(is_mutated=None)
        fragment_list.append(unknown_fragment)

        self.assertEqual(fragment_list.count_mutated_fragments(), 1)

    def test_fragment_list_count_endmotifs(self) -> None:
        """Test end motif counting."""
        fragment_list = FragmentList()

        fragment1 = create_mock_fragment(
            end5p="ATCG", end3p="CGAT", is_bogus=False
        )
        fragment2 = create_mock_fragment(
            end5p="ATCG", end3p="TGCA", is_bogus=False
        )
        fragment3 = create_mock_fragment(
            end5p="NNNN", end3p="CGAT", is_bogus=True
        )

        fragment_list.append(fragment1)
        fragment_list.append(fragment2)
        fragment_list.append(fragment3)
        motifs_5p, motifs_3p, num_frags = fragment_list.count_endmotifs(3)

        self.assertEqual(num_frags, 2)
        self.assertEqual(motifs_5p["ATC"], 2)
        self.assertEqual(motifs_3p["GAT"], 1)
        self.assertEqual(motifs_3p["GCA"], 1)

    def test_fragment_list_to_frag_file(self) -> None:
        """Test saving FragmentList to .frag file."""
        fragment_list = create_test_fragment_list(5)

        with tempfile.TemporaryDirectory() as temp_dir:
            fragment_list.to_frag_file("test", temp_dir)
            frag_path = os.path.join(temp_dir, "test.frag")
            self.assertTrue(os.path.exists(frag_path))

    def test_fragment_list_to_interval_table(self) -> None:
        """Test conversion to IntervalTable."""
        fragment_list = FragmentList()

        frag1 = create_mock_fragment(
            start_pos=100, end_pos=200, chrom="1", is_bogus=False
        )
        frag2 = create_mock_fragment(
            start_pos=300, end_pos=400, chrom="1", is_bogus=False
        )
        frag3 = create_mock_fragment(
            start_pos=500, end_pos=600, chrom="2", is_bogus=False
        )
        frag4 = create_mock_fragment(
            start_pos=700, end_pos=800, chrom="1", is_bogus=True
        )
        fragment_list.append(frag1)
        fragment_list.append(frag2)
        fragment_list.append(frag3)
        fragment_list.append(frag4)
        interval_table = fragment_list.to_interval_table()

        overlaps_chr1: list[Fragment] = \
            interval_table.get_overlaps("1", 150, 350)
        self.assertEqual(len(overlaps_chr1), 2)

        overlaps_chr2: list[Fragment] = \
            interval_table.get_overlaps("2", 550, 650)
        self.assertEqual(len(overlaps_chr2), 1)


class TestFragmentCollection(unittest.TestCase):
    """Test FragmentCollection class functionality."""

    def test_fragment_collection_basic_operations(self) -> None:
        """Test basic FragmentCollection operations."""
        collection = FragmentCollection()

        fragment_list1 = create_test_fragment_list(5)
        fragment_list2 = create_test_fragment_list(3)
        collection.append("sample1", fragment_list1)
        collection.append("sample2", fragment_list2)

        names = collection.record_names()
        self.assertIn("sample1", names)
        self.assertIn("sample2", names)

        items = list(collection)
        self.assertEqual(len(items), 2)

    def test_fragment_collection_to_frag_files(self) -> None:
        """Test saving FragmentCollection to .frag files."""
        collection = FragmentCollection()

        fragment_list1 = create_test_fragment_list(3)
        fragment_list2 = create_test_fragment_list(2)

        collection.append("sample1", fragment_list1)
        collection.append("sample2", fragment_list2)

        with tempfile.TemporaryDirectory() as temp_dir:
            collection.to_frag_files(temp_dir)
            self.assertTrue(
                os.path.exists(os.path.join(temp_dir, "sample1.frag"))
            )
            self.assertTrue(
                os.path.exists(os.path.join(temp_dir, "sample2.frag"))
            )


class TestIntervalTable(unittest.TestCase):
    """Test IntervalTable class functionality."""

    def test_interval_table_insert_and_query(self) -> None:
        """Test IntervalTable insertion and querying."""
        table = IntervalTable()

        frag1 = create_mock_fragment(start_pos=100, end_pos=200, chrom="1")
        frag2 = create_mock_fragment(start_pos=300, end_pos=400, chrom="1")
        frag3 = create_mock_fragment(start_pos=150, end_pos=250, chrom="2")
        table.insert(frag1)
        table.insert(frag2)
        table.insert(frag3)

        overlaps: list[Fragment] = table.get_overlaps("1", 150, 350)
        self.assertEqual(len(overlaps), 2)

        overlaps = table.get_overlaps("2", 200, 300)
        self.assertEqual(len(overlaps), 1)

        overlaps = table.get_overlaps("1", 500, 600)
        self.assertEqual(len(overlaps), 0)


class TestUtilityFunctions(unittest.TestCase):
    """Test utility functions in fragment module."""

    def test_is_duplex(self) -> None:
        """Test duplex read detection."""
        duplex_read = Mock()

        def duplex_get_tag(tag: str) -> str:
            if tag in ["XI", "XJ"]:
                return "value"
            raise KeyError(f"Tag {tag} not found")

        duplex_read.get_tag.side_effect = duplex_get_tag  # type: ignore
        self.assertTrue(is_duplex(duplex_read))

        non_duplex_read = Mock()
        non_duplex_read.get_tag.side_effect = (  # type: ignore
            KeyError("No tags found")
        )
        self.assertFalse(is_duplex(non_duplex_read))

        partial_duplex_read = Mock()

        def partial_get_tag(tag: str) -> str:
            if tag == "XI":
                return "value"
            raise KeyError(f"Tag {tag} not found")

        partial_duplex_read.get_tag.side_effect = (  # type: ignore
            partial_get_tag
        )

        self.assertFalse(is_duplex(partial_duplex_read))


class TestVCFParsing(unittest.TestCase):
    """Test VCF parsing and mutation annotation functionality."""

    def test_build_mutated_reads_set_basic(self) -> None:
        """Test basic VCF parsing with simple SNV."""
        bam_file = Mock()
        bam_file.header.to_dict.return_value = {  # type: ignore
            "SQ": [{"SN": "chr1"}]  # type: ignore
        }

        vcf_file = Mock()
        vcf_file.filename = b"test.vcf"
        variant = create_mock_variant_record(
            contig="chr1", start=125, stop=126, ref="A", alts=("T",)
        )
        vcf_file.__iter__ = Mock(return_value=iter([variant]))

        mutated_read = MockAlignedSegment(
            query_name="read1",
            query_sequence="ATCGATCGTCG",  # T at position 5 (alt)
            reference_start=120,
            reference_end=131
        )
        mutated_read.get_reference_positions = Mock(  # type: ignore
            return_value=list(range(120, 131))  # type: ignore
        )

        normal_read = MockAlignedSegment(
            query_name="read2",
            query_sequence="ATCGAACGTCG",  # A at position 5 (ref)
            reference_start=120,
            reference_end=131
        )
        normal_read.get_reference_positions = Mock(  # type: ignore
            return_value=list(range(120, 131))  # type: ignore
        )

        bam_file.fetch.return_value = [  # type: ignore
            mutated_read, normal_read
        ]

        result = Fragment.build_mutated_reads_set(bam_file, vcf_file)

        self.assertEqual(len(result), 1)
        self.assertIn(mutated_read.query_name, result)
        self.assertNotIn(normal_read.query_name, result)

    def test_build_mutated_reads_set_non_snv_filtered(self) -> None:
        """Test that non-SNV variants are filtered out."""
        bam_file = Mock()
        bam_file.header.to_dict.return_value = {  # type: ignore
            "SQ": [{"SN": "chr1"}]  # type: ignore
        }

        vcf_file = Mock()
        vcf_file.filename = b"test.vcf"
        variant = create_mock_variant_record(
            contig="chr1", start=125, stop=127, ref="AT", alts=("A",), rlen=2
        )
        vcf_file.__iter__ = Mock(return_value=iter([variant]))

        result = Fragment.build_mutated_reads_set(bam_file, vcf_file)

        self.assertEqual(len(result), 0)
        bam_file.fetch.assert_not_called()  # type: ignore

    def test_build_mutated_reads_set_chromosome_name_normalization(
        self
    ) -> None:
        """Test chromosome name normalization between VCF and BAM."""
        bam_file = Mock()
        bam_file.header.to_dict.return_value = {  # type: ignore
            "SQ": [{"SN": "chr1"}]  # type: ignore
        }

        vcf_file = Mock()
        vcf_file.filename = b"test.vcf"
        variant = create_mock_variant_record(
            contig="1", start=125, stop=126, ref="A", alts=("T",)
        )
        vcf_file.__iter__ = Mock(return_value=iter([variant]))

        mutated_read = MockAlignedSegment(
            query_name="read1",
            query_sequence="ATCGATCGTCG",  # T at position 5 (alt)
            reference_start=120,
            reference_end=131
        )
        mutated_read.get_reference_positions = Mock(  # type: ignore
            return_value=list(range(120, 131))  # type: ignore
        )
        bam_file.fetch.return_value = [mutated_read]  # type: ignore

        result = Fragment.build_mutated_reads_set(bam_file, vcf_file)
        bam_file.fetch.assert_called_with(  # type: ignore
            contig="chr1", start=125, stop=126
        )
        self.assertEqual(len(result), 1)

    def test_build_mutated_reads_set_read_position_not_in_variant(
        self
    ) -> None:
        """Test reads that do not span the variant position."""
        bam_file = Mock()
        bam_file.header.to_dict.return_value = {  # type: ignore
            "SQ": [{"SN": "chr1"}]  # type: ignore
        }

        vcf_file = Mock()
        vcf_file.filename = b"test.vcf"
        variant = create_mock_variant_record(
            contig="chr1", start=125, stop=126, ref="A", alts=("T",)
        )
        vcf_file.__iter__ = Mock(return_value=iter([variant]))

        read = MockAlignedSegment(
            query_name="read1",
            query_sequence="ATCGTACGTCG",
            reference_start=130,
            reference_end=141
        )
        read.get_reference_positions = Mock(  # type: ignore
            return_value=list(range(130, 141))  # type: ignore
        )

        bam_file.fetch.return_value = [read]  # type: ignore
        result = Fragment.build_mutated_reads_set(bam_file, vcf_file)

        self.assertEqual(len(result), 0)

    def test_build_mutated_reads_set_unknown_variant_base(self) -> None:
        """Test reads with bases that are neither ref nor alt."""
        bam_file = Mock()
        bam_file.header.to_dict.return_value = {  # type: ignore
            "SQ": [{"SN": "chr1"}]  # type: ignore
        }

        vcf_file = Mock()
        vcf_file.filename = b"test.vcf"
        variant = create_mock_variant_record(
            contig="chr1", start=125, stop=126, ref="A", alts=("T",)
        )
        vcf_file.__iter__ = Mock(return_value=iter([variant]))

        read = MockAlignedSegment(
            query_name="read1",
            query_sequence="ATCGGACGTCG",  # G at position 5
            reference_start=120,
            reference_end=131
        )
        read.get_reference_positions = Mock(  # type: ignore
            return_value=list(range(120, 131))  # type: ignore
        )

        bam_file.fetch.return_value = [read]  # type: ignore
        result = Fragment.build_mutated_reads_set(bam_file, vcf_file)

        self.assertEqual(len(result), 0)

    def test_build_mutated_reads_set_multiple_variants(self) -> None:
        """Test processing multiple variants in VCF."""
        bam_file = Mock()
        bam_file.header.to_dict.return_value = {  # type: ignore
            "SQ": [{"SN": "chr1"}]  # type: ignore
        }

        vcf_file = Mock()
        vcf_file.filename = b"test.vcf"
        variant1 = create_mock_variant_record(
            contig="chr1", start=125, stop=126, ref="A", alts=("T",)
        )
        variant2 = create_mock_variant_record(
            contig="chr1", start=225, stop=226, ref="C", alts=("G",)
        )
        vcf_file.__iter__ = Mock(return_value=iter([variant1, variant2]))

        read1 = MockAlignedSegment(
            query_name="read1", query_sequence="ATCGATCGTCG"  # T at position 5
        )
        read1.get_reference_positions = Mock(  # type: ignore
            return_value=list(range(120, 131))  # type: ignore
        )

        read2 = MockAlignedSegment(
            query_name="read2", query_sequence="ATCGGGCGTCG"  # G at position 5
        )
        read2.get_reference_positions = Mock(  # type: ignore
            return_value=list(range(220, 231))  # type: ignore
        )

        bam_file.fetch.side_effect = [[read1], [read2]]  # type: ignore
        result = Fragment.build_mutated_reads_set(bam_file, vcf_file)

        self.assertEqual(len(result), 2)
        self.assertIn(read1.query_name, result)
        self.assertIn(read2.query_name, result)

    def test_build_mutated_reads_set_invalid_chromosome_error(self) -> None:
        """Test error handling for invalid chromosome in VCF."""
        bam_file = Mock()
        bam_file.header.to_dict.return_value = {  # type: ignore
            "SQ": [{"SN": "chr1"}]  # type: ignore
        }

        vcf_file = Mock()
        vcf_file.filename = b"test.vcf"
        variant = create_mock_variant_record(
            contig="invalid_chr", start=125, stop=126, ref="A", alts=("T",)
        )
        vcf_file.__iter__ = Mock(return_value=iter([variant]))

        bam_file.fetch.side_effect = ValueError(  # type: ignore
            "Invalid chromosome"
        )

        with self.assertRaises(SystemExit):
            Fragment.build_mutated_reads_set(bam_file, vcf_file)


class TestVCFIntegration(unittest.TestCase):
    """Integration tests for VCF parsing with real BAM processing."""

    def test_fragment_mutation_consistency(self) -> None:
        """Test that mutation annotation is consistent across fragment ops."""
        fragments = FragmentList()

        mutated_frag = create_mock_fragment(is_mutated=True)
        normal_frag = create_mock_fragment(is_mutated=False)
        unknown_frag = create_mock_fragment(is_mutated=None)

        fragments.append(mutated_frag)
        fragments.append(normal_frag)
        fragments.append(unknown_frag)

        self.assertEqual(fragments.count_mutated_fragments(), 1)
        self.assertEqual(fragments.length(), 3)

        temp_path = create_temp_frag_file(fragments)
        try:
            from pyfraglib.fragfile import FragFile
            frag_file = FragFile(temp_path)
            loaded_fragments = frag_file.get_fragment_list()

            self.assertEqual(loaded_fragments.count_mutated_fragments(), 1)
            for orig, loaded in zip(fragments, loaded_fragments):
                self.assertEqual(orig.is_mutated, loaded.is_mutated)

            frag_file.close()
        finally:
            cleanup_temp_file(temp_path)

    @unittest.skipUnless(  # type: ignore
        os.path.exists(MUTATED_READS_ONLY_TEST_BAM) and
        os.path.exists(MUTATED_READS_ONLY_TEST_BAM_BAI) and
        os.path.exists(MUTATED_READS_ONLY_TEST_VCF),
        "Integration test files not found"
    )
    def test_mutated_only_bam_integration(self) -> None:
        """
        Test that fragments in mutated-only BAM are correctly annotated. This
        integration test relies on sensitive data and is thus ignored in most
        contexts.
        """
        fragments = Fragment.from_bam(
            MUTATED_READS_ONLY_TEST_BAM, MUTATED_READS_ONLY_TEST_VCF
        )
        self.assertGreater(
            fragments.length(), 0, "No fragments loaded from BAM file"
        )

        total_fragments: int = fragments.length()
        mutated_fragments: int = fragments.count_mutated_fragments()
        mutated_count: int = 0
        non_mutated_count: int = 0
        unknown_count: int = 0

        for fragment in fragments:
            if fragment.is_mutated is True:
                mutated_count += 1
            elif fragment.is_mutated is False:
                non_mutated_count += 1
            else:  # fragment.is_mutated is None
                unknown_count += 1

        self.assertEqual(
            mutated_count, mutated_fragments,
            "count_mutated_fragments() should match manual count"
        )

        mutation_rate = mutated_fragments / total_fragments
        non_mutation_rate = non_mutated_count / total_fragments
        known_fragments = mutated_count + non_mutated_count
        known_rate = known_fragments / total_fragments
        self.assertEqual(
            mutation_rate, 1.0,
            f"Expected all fragments to be mutated, "
            f"but only {mutation_rate:.2%} "
            f"({mutated_fragments}/{total_fragments}) were mutated"
        )
        self.assertEqual(
            non_mutation_rate, 0.0,
            f"Expected no fragments to be non-mutated, "
            f"but {non_mutation_rate:.2%} "
            f"({non_mutated_count}/{total_fragments}) were non-mutated"
        )
        self.assertEqual(
            known_rate, 1.0,
            f"Expected all fragments to have known "
            f"mutation status, but only {known_rate:.2%} "
            f"({known_fragments}/{total_fragments}) were known"
        )

    @unittest.skipUnless(  # type: ignore
        os.path.exists(UNMUTATED_READS_ONLY_TEST_BAM) and
        os.path.exists(UNMUTATED_READS_ONLY_TEST_BAM_BAI) and
        os.path.exists(UNMUTATED_READS_ONLY_TEST_VCF),
        "Integration test files not found"
    )
    def test_unmutated_only_bam_integration(self) -> None:
        """
        Test that fragments in unmutated-only BAM are correctly annotated.
        This integration test complements the mutated-only test and relies
        on sensitive data, thus is ignored in most contexts.
        """
        fragments = Fragment.from_bam(
            UNMUTATED_READS_ONLY_TEST_BAM, UNMUTATED_READS_ONLY_TEST_VCF
        )
        self.assertGreater(
            fragments.length(), 0, "No fragments loaded from BAM file"
        )

        total_fragments: int = fragments.length()
        mutated_fragments: int = fragments.count_mutated_fragments()
        mutated_count: int = 0
        non_mutated_count: int = 0
        unknown_count: int = 0

        for fragment in fragments:
            if fragment.is_mutated is True:
                mutated_count += 1
            elif fragment.is_mutated is False:
                non_mutated_count += 1
            else:  # fragment.is_mutated is None
                unknown_count += 1

        self.assertEqual(
            mutated_count, mutated_fragments,
            "count_mutated_fragments() should match manual count"
        )

        mutation_rate = mutated_fragments / total_fragments
        non_mutation_rate = non_mutated_count / total_fragments
        known_fragments = mutated_count + non_mutated_count
        known_rate = known_fragments / total_fragments

        self.assertEqual(
            mutation_rate, 0.0,
            f"Expected no fragments to be mutated, "
            f"but {mutation_rate:.2%} "
            f"({mutated_fragments}/{total_fragments}) were mutated"
        )
        self.assertEqual(
            non_mutation_rate, 1.0,
            f"Expected all fragments to be non-mutated, "
            f"but only {non_mutation_rate:.2%} "
            f"({non_mutated_count}/{total_fragments}) were non-mutated"
        )
        self.assertEqual(
            known_rate, 1.0,
            f"Expected all fragments to have known "
            f"mutation status, but only {known_rate:.2%} "
            f"({known_fragments}/{total_fragments}) were known"
        )

    @unittest.skipUnless(  # type: ignore
        os.path.exists(MIXED_MUTATED_READS_TEST_BAM) and
        os.path.exists(MIXED_MUTATED_READS_TEST_BAM_BAI) and
        os.path.exists(MIXED_MUTATED_READS_TEST_VCF),
        "Integration test files not found"
    )
    def test_mixed_mutated_bam_integration(self) -> None:
        """
        Test that fragments in BAM with mutationally mixed reads are correctly
        annotated. Compare to the mutated-/unmutated-reads-only tests for more
        information.
        """
        fragments = Fragment.from_bam(
            MIXED_MUTATED_READS_TEST_BAM, MIXED_MUTATED_READS_TEST_VCF
        )
        self.assertGreater(
            fragments.length(), 0, "No fragments loaded from BAM file"
        )

        total_fragments: int = fragments.length()
        mutated_fragments: int = fragments.count_mutated_fragments()
        mutated_count: int = 0
        non_mutated_count: int = 0
        unknown_count: int = 0

        for fragment in fragments:
            if fragment.is_mutated is True:
                mutated_count += 1
            elif fragment.is_mutated is False:
                non_mutated_count += 1
            else:  # fragment.is_mutated is None
                unknown_count += 1

        self.assertEqual(
            mutated_count, mutated_fragments,
            "count_mutated_fragments() should match manual count"
        )

        mutation_rate = mutated_fragments / total_fragments
        non_mutation_rate = non_mutated_count / total_fragments
        known_fragments = mutated_count + non_mutated_count
        known_rate = known_fragments / total_fragments

        self.assertNotEqual(
            mutation_rate, 0.0,
            f"Expected some fragments to be mutated, "
            f"but none were ({mutation_rate:.2%}, "
            f"{mutated_fragments}/{total_fragments})"
        )
        self.assertNotEqual(
            non_mutation_rate, 1.0,
            f"Expected not all fragments to be non-mutated, "
            f"but all were ({non_mutation_rate:.2%}, "
            f"({non_mutated_count}/{total_fragments})"
        )
        self.assertEqual(
            known_rate, 1.0,
            f"Expected all fragments to have known "
            f"mutation status, but only {known_rate:.2%} "
            f"({known_fragments}/{total_fragments}) were known"
        )


class TestVCFEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions in VCF processing."""

    def test_empty_vcf_file(self) -> None:
        """Test handling of empty VCF file."""
        bam_file = Mock()
        vcf_file = Mock()
        vcf_file.filename = b"empty.vcf"
        vcf_file.__iter__ = Mock(return_value=iter([]))

        result = Fragment.build_mutated_reads_set(bam_file, vcf_file)

        self.assertEqual(len(result), 0)

    def test_vcf_with_multiple_alt_alleles(self) -> None:
        """Test VCF records with multiple alternative alleles."""
        bam_file = Mock()
        bam_file.header.to_dict.return_value = {  # type: ignore
            "SQ": [{"SN": "chr1"}]  # type: ignore
        }

        vcf_file = Mock()
        vcf_file.filename = b"test.vcf"
        variant = create_mock_variant_record(
            contig="chr1", start=125, stop=126, ref="A", alts=("T", "G")
        )
        vcf_file.__iter__ = Mock(return_value=iter([variant]))

        read1 = MockAlignedSegment(
            query_name="read1",
            query_sequence="ATCGATCGTCG",  # T at position 5 (alt)
            reference_start=120,
            reference_end=131
        )
        read1.get_reference_positions = Mock(  # type: ignore
            return_value=list(range(120, 131))  # type: ignore
        )

        read2 = MockAlignedSegment(
            query_name="read2",
            query_sequence="ATCGGGCGTCG",  # G at position 5 (alt)
            reference_start=120,
            reference_end=131
        )
        read2.get_reference_positions = Mock(  # type: ignore
            return_value=list(range(120, 131))  # type: ignore
        )

        bam_file.fetch.return_value = [read1, read2]  # type: ignore
        result = Fragment.build_mutated_reads_set(bam_file, vcf_file)

        self.assertEqual(len(result), 2)
        self.assertIn(read1.query_name, result)
        self.assertIn(read2.query_name, result)

    def test_read_with_none_query_sequence(self) -> None:
        """Test handling of reads with None query sequence."""
        bam_file = Mock()
        bam_file.header.to_dict.return_value = {  # type: ignore
            "SQ": [{"SN": "chr1"}]  # type: ignore
        }

        vcf_file = Mock()
        vcf_file.filename = b"test.vcf"
        variant = create_mock_variant_record(
            contig="chr1", start=125, stop=126, ref="A", alts=("T",)
        )
        vcf_file.__iter__ = Mock(return_value=iter([variant]))

        read = MockAlignedSegment(
            query_name="read1",
            query_sequence=None,
            reference_start=120,
            reference_end=131
        )
        read.get_reference_positions = Mock(  # type: ignore
            return_value=list(range(120, 131))  # type: ignore
        )
        bam_file.fetch.return_value = [read]  # type: ignore

        with self.assertRaises(AssertionError):
            Fragment.build_mutated_reads_set(bam_file, vcf_file)


class TestConstants(unittest.TestCase):
    """Test module constants."""

    def test_valid_chromosome_names(self) -> None:
        """Test valid chromosome names list."""
        self.assertIn("1", VALID_CHROMOSOME_NAMES)
        self.assertIn("22", VALID_CHROMOSOME_NAMES)
        self.assertIn("X", VALID_CHROMOSOME_NAMES)
        self.assertIn("Y", VALID_CHROMOSOME_NAMES)
        self.assertIn("chr1", VALID_CHROMOSOME_NAMES)
        self.assertIn("chrX", VALID_CHROMOSOME_NAMES)
        self.assertIn("chrY", VALID_CHROMOSOME_NAMES)
        self.assertIn("M", VALID_CHROMOSOME_NAMES)
        self.assertIn("m", VALID_CHROMOSOME_NAMES)

    def test_constants_values(self) -> None:
        """Test that constants have expected values."""
        self.assertEqual(INSERT_SIZE_UPPER_BOUND, 900)
        self.assertEqual(MIN_MAPQ, 20)
        self.assertEqual(DEFAULT_KMER_LEN, 3)
        self.assertEqual(MAX_KMER_LEN, 4)


if __name__ == "__main__":
    unittest.main()
