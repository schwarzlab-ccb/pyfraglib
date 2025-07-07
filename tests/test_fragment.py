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
from pyfraglib.fragment import Fragment, FragmentList, FragmentCollection, \
                               IntervalTable, is_duplex, \
                               VALID_CHROMOSOME_NAMES, \
                               INSERT_SIZE_UPPER_BOUND, MIN_MAPQ, \
                               DEFAULT_KMER_LEN, MAX_KMER_LEN
from tests.test_fixtures import MockAlignedSegment, create_mock_fragment, \
                                create_test_fragment_list


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
        mutated_reads: set[MockAlignedSegment] = {read}

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
