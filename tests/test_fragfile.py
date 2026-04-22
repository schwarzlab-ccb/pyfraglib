# This file is part of `pyfraglib`, a software suite to calculate fragmentomics
# features from cfDNA and perform downstream analyses.
#
# Copyright (C) 2026 Daniel Schütte, daniel.schuette@iccb-cologne.org
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

from pyfraglib.fragfile import FragFile
from pyfraglib.fragment import Fragment, FragmentList
from tests.test_fixtures import create_mock_fragment, \
                                create_test_fragment_list, \
                                create_temp_frag_file, cleanup_temp_file


class TestFragFile(unittest.TestCase):
    """Test FragFile class functionality."""

    def setUp(self) -> None:
        """Set up test fixtures."""
        self.test_fragments = create_test_fragment_list(5)
        self.temp_frag_path = create_temp_frag_file(self.test_fragments)

    def tearDown(self) -> None:
        """Clean up test files."""
        cleanup_temp_file(self.temp_frag_path)

    def test_fragfile_init(self) -> None:
        """Test FragFile initialization."""
        frag_file = FragFile(self.temp_frag_path)
        self.assertEqual(
            frag_file._FragFile__path, self.temp_frag_path  # type: ignore
        )
        self.assertFalse(frag_file.closed)
        frag_file.close()

    def test_fragfile_iteration(self) -> None:
        """Test FragFile iteration over fragments."""
        frag_file = FragFile(self.temp_frag_path)

        fragments: list[Fragment] = []
        for fragment in frag_file:
            fragments.append(fragment)

        self.assertEqual(len(fragments), 5)

        for i, fragment in enumerate(fragments):
            self.assertIsInstance(fragment, Fragment)
            self.assertEqual(fragment.start_pos, 100 + i * 200)
            self.assertEqual(fragment.end_pos, 250 + i * 200)
            self.assertEqual(fragment.length, 150 + i * 10)

        frag_file.close()

    def test_fragfile_get_fragment_list(self) -> None:
        """Test getting complete FragmentList from FragFile."""
        frag_file = FragFile(self.temp_frag_path)

        fragment_list = frag_file.get_fragment_list()

        self.assertIsInstance(fragment_list, FragmentList)
        self.assertEqual(fragment_list.length(), 5)

        fragments = list(fragment_list)
        for i, fragment in enumerate(fragments):
            self.assertEqual(fragment.start_pos, 100 + i * 200)
            self.assertEqual(fragment.end_pos, 250 + i * 200)
            self.assertEqual(fragment.length, 150 + i * 10)

        frag_file.close()

    def test_fragfile_close(self) -> None:
        """Test FragFile close method."""
        frag_file = FragFile(self.temp_frag_path)
        self.assertFalse(frag_file.closed)

        frag_file.close()
        self.assertTrue(frag_file.closed)

    def test_fragfile_read_after_close_raises(self) -> None:
        """Reads after close must raise rather than returning silently."""
        frag_file = FragFile(self.temp_frag_path)
        frag_file.close()
        with self.assertRaises(ValueError):
            list(frag_file)
        with self.assertRaises(ValueError):
            frag_file.get_fragment_list()

    def test_fragfile_reopen_after_close(self) -> None:
        """After closing, a fresh FragFile should read the same data."""
        frag_file = FragFile(self.temp_frag_path)
        fragments1 = list(frag_file)
        self.assertEqual(len(fragments1), 5)
        frag_file.close()

        frag_file2 = FragFile(self.temp_frag_path)
        fragments2 = list(frag_file2)
        self.assertEqual(len(fragments2), 5)
        frag_file2.close()

    def test_fragfile_empty_file(self) -> None:
        """Test FragFile with an empty fragment list."""
        empty_list = FragmentList()
        empty_path = create_temp_frag_file(empty_list)

        try:
            frag_file = FragFile(empty_path)
            fragments = list(frag_file)
            self.assertEqual(len(fragments), 0)
            frag_file.close()

            frag_file = FragFile(empty_path)
            fragment_list = frag_file.get_fragment_list()
            self.assertEqual(fragment_list.length(), 0)
            frag_file.close()
        finally:
            cleanup_temp_file(empty_path)

    def test_fragfile_single_fragment(self) -> None:
        """Test FragFile with a single fragment."""
        fragments = FragmentList()
        fragments.append(create_mock_fragment(
            start_pos=500, end_pos=650, chrom="X", length=150,
        ))
        single_frag_path = create_temp_frag_file(fragments)

        try:
            frag_file = FragFile(single_frag_path)
            loaded = list(frag_file)
            self.assertEqual(len(loaded), 1)

            fragment = loaded[0]
            self.assertEqual(fragment.start_pos, 500)
            self.assertEqual(fragment.end_pos, 650)
            self.assertEqual(fragment.chrom, "X")
            self.assertEqual(fragment.length, 150)

            frag_file.close()
        finally:
            cleanup_temp_file(single_frag_path)

    def test_fragfile_nonexistent_file(self) -> None:
        """Test FragFile with a nonexistent file."""
        with self.assertRaises(FileNotFoundError):
            FragFile("/nonexistent/path/test.frag")

    def test_fragfile_corrupted_file(self) -> None:
        """Test FragFile with a corrupted (non-Parquet) .frag file."""
        with tempfile.NamedTemporaryFile(suffix=".frag", delete=False) as tmp:
            tmp.write(b"This is not a valid Parquet file")
            corrupted_path = tmp.name

        try:
            with self.assertRaises(Exception):
                FragFile(corrupted_path)
        finally:
            cleanup_temp_file(corrupted_path)

    def test_fragfile_destructor(self) -> None:
        """Test FragFile destructor closes file."""
        frag_file = FragFile(self.temp_frag_path)
        self.assertFalse(frag_file.closed)
        del frag_file
        self.assertTrue(True)  # destructor did not raise


class TestFragFileIntegration(unittest.TestCase):
    """Test FragFile integration with other components."""

    def test_fragfile_with_varied_fragments(self) -> None:
        """Test FragFile with fragments having varied properties."""
        fragments = FragmentList()
        normal_frag = create_mock_fragment(
            start_pos=100, end_pos=250, chrom="1",
            is_bogus=False, is_mutated=False
        )
        fragments.append(normal_frag)

        bogus_frag = create_mock_fragment(
            start_pos=300, end_pos=450, chrom="2",
            is_bogus=True, is_mutated=None
        )
        fragments.append(bogus_frag)

        mutated_frag = create_mock_fragment(
            start_pos=500, end_pos=650, chrom="X",
            is_bogus=False, is_mutated=True
        )
        fragments.append(mutated_frag)

        temp_path = create_temp_frag_file(fragments)
        try:
            frag_file = FragFile(temp_path)
            loaded_fragments = frag_file.get_fragment_list()

            self.assertEqual(loaded_fragments.length(), 3)
            self.assertEqual(loaded_fragments.count_bogus_fragments(), 1)
            self.assertEqual(loaded_fragments.count_mutated_fragments(), 1)

            frag_file.close()
        finally:
            cleanup_temp_file(temp_path)

    def test_fragfile_roundtrip_consistency(self) -> None:
        """Test that fragments survive a save/load cycle unchanged."""
        original_fragments = create_test_fragment_list(10)
        temp_path = create_temp_frag_file(original_fragments)

        try:
            frag_file = FragFile(temp_path)
            loaded_fragments = frag_file.get_fragment_list()
            self.assertEqual(
                original_fragments.length(), loaded_fragments.length()
            )

            original_list = list(original_fragments)
            loaded_list = list(loaded_fragments)

            for orig, loaded in zip(original_list, loaded_list):
                self.assertEqual(orig.start_pos, loaded.start_pos)
                self.assertEqual(orig.end_pos, loaded.end_pos)
                self.assertEqual(orig.chrom, loaded.chrom)
                self.assertEqual(orig.length, loaded.length)
                self.assertEqual(orig.end5p, loaded.end5p)
                self.assertEqual(orig.end3p, loaded.end3p)
                self.assertEqual(orig.is_bogus, loaded.is_bogus)
                self.assertEqual(orig.is_mutated, loaded.is_mutated)
                self.assertEqual(orig.is_forward, loaded.is_forward)

            frag_file.close()
        finally:
            cleanup_temp_file(temp_path)


if __name__ == "__main__":
    unittest.main()
