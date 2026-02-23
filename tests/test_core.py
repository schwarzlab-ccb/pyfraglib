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
import os

from unittest.mock import patch
from pyfraglib.core import get_chromosome_length, homogenize_contig_name, \
                           homogenize_to_chrom_naming_convention, \
                           shannon_entropy, simpson_index, detect_cpus, \
                           hg19_chromosomes, hg38_chromosomes, \
                           PyfraglibException, CodeUnreachableError


class TestChromosomeFunctions(unittest.TestCase):
    """Test chromosome-related functions."""

    def test_get_chromosome_length_hg19(self) -> None:
        """Test getting chromosome length for hg19."""
        length = get_chromosome_length("1", "hg19")
        self.assertEqual(length, 249250621)

        length = get_chromosome_length("chr1", "hg19")
        self.assertEqual(length, 249250621)

        length_x = get_chromosome_length("X", "hg19")
        self.assertEqual(length_x, 155270560)

        length_y = get_chromosome_length("Y", "hg19")
        self.assertEqual(length_y, 59373566)

        length_m = get_chromosome_length("M", "hg19")
        self.assertEqual(length_m, 16569)

    def test_get_chromosome_length_hg38(self) -> None:
        """Test getting chromosome length for hg38."""
        length = get_chromosome_length("1", "hg38")
        self.assertEqual(length, 248956422)

        length_hg19 = get_chromosome_length("1", "hg19")
        self.assertNotEqual(length, length_hg19)

    def test_get_chromosome_length_invalid_genome(self) -> None:
        """Test that invalid genome raises SystemExit."""
        with self.assertRaises(SystemExit):
            get_chromosome_length("1", "invalid_genome")

    def test_get_chromosome_length_invalid_chromosome(self) -> None:
        """Test that invalid chromosome raises SystemExit."""
        with self.assertRaises(SystemExit):
            get_chromosome_length("invalid_chr", "hg19")

    def test_homogenize_contig_name(self) -> None:
        """Test chromosome name homogenization."""
        self.assertEqual(homogenize_contig_name("chr1"), "1")
        self.assertEqual(homogenize_contig_name("chrX"), "X")
        self.assertEqual(homogenize_contig_name("chrY"), "Y")
        self.assertEqual(homogenize_contig_name("1"), "1")
        self.assertEqual(homogenize_contig_name("X"), "X")
        self.assertEqual(homogenize_contig_name("M"), "M")

    def test_homogenize_to_chrom_naming_convention(self) -> None:
        """Test chromosome naming convention homogenization."""
        header_chr: dict[str, object] = {"SQ": [{"SN": "chr1"}]}
        result = homogenize_to_chrom_naming_convention("1", header_chr)
        self.assertEqual(result, "chr1")

        result = homogenize_to_chrom_naming_convention("chr1", header_chr)
        self.assertEqual(result, "chr1")

        header_no_chr: dict[str, object] = {"SQ": [{"SN": "1"}]}
        result = homogenize_to_chrom_naming_convention("chr1", header_no_chr)
        self.assertEqual(result, "1")

        result = homogenize_to_chrom_naming_convention("1", header_no_chr)
        self.assertEqual(result, "1")

        invalid_header: dict[str, object] = {}
        result = homogenize_to_chrom_naming_convention("chr1", invalid_header)
        self.assertEqual(result, "chr1")


class TestDiversityIndices(unittest.TestCase):
    """Test diversity index calculations."""

    def test_shannon_entropy(self) -> None:
        """Test Shannon entropy calculation."""
        uniform_props = [0.25, 0.25, 0.25, 0.25]
        entropy = shannon_entropy(uniform_props)
        self.assertAlmostEqual(entropy, 1.386, places=3)

        dominant_props = [0.97, 0.01, 0.01, 0.01]
        entropy = shannon_entropy(dominant_props)
        self.assertLess(entropy, 0.2)

        entropy = shannon_entropy([])
        self.assertEqual(entropy, 0.0)

    def test_simpson_index(self) -> None:
        """Test Simpson index calculation."""
        uniform_props = [0.25, 0.25, 0.25, 0.25]
        simpson = simpson_index(uniform_props)
        self.assertAlmostEqual(simpson, 0.25, places=3)

        dominant_props = [0.97, 0.01, 0.01, 0.01]
        simpson = simpson_index(dominant_props)
        self.assertGreater(simpson, 0.9)

        simpson = simpson_index([])
        self.assertEqual(simpson, 0.0)


class TestUtilityFunctions(unittest.TestCase):
    """Test utility functions."""

    @patch.dict(os.environ, {"SLURM_CPUS_PER_TASK": "8"})
    def test_detect_cpus_with_slurm(self) -> None:
        """Test CPU detection with SLURM environment variable."""
        cpus = detect_cpus()
        self.assertEqual(cpus, 8)

    @patch.dict(os.environ, {}, clear=True)
    def test_detect_cpus_without_slurm(self) -> None:
        """Test CPU detection without SLURM environment variable."""
        cpus = detect_cpus()
        cpus_without_slurm = os.cpu_count() if os.cpu_count is not None else 1
        self.assertEqual(cpus, cpus_without_slurm)


class TestExceptions(unittest.TestCase):
    """Test custom exceptions."""

    def test_pyfraglib_exception(self) -> None:
        """Test PyfraglibException creation."""
        with self.assertLogs("pyfraglib", level="FATAL"):
            exception = PyfraglibException("test error message")
            self.assertIsInstance(exception, Exception)

    def test_code_unreachable_error(self) -> None:
        """Test CodeUnreachableError creation."""
        with self.assertLogs("pyfraglib", level="FATAL"):
            error = CodeUnreachableError("unreachable code")
            self.assertIsInstance(error, PyfraglibException)


class TestChromosomeData(unittest.TestCase):
    """Test chromosome data constants."""

    def test_hg19_chromosomes_structure(self) -> None:
        """Test hg19 chromosome data structure."""
        self.assertEqual(len(hg19_chromosomes), 25)  # 22 autosomes + X, Y, M

        chr1 = hg19_chromosomes[0]
        self.assertEqual(chr1[0], "1")  # chromosome name
        self.assertEqual(chr1[1], 249250621)  # length
        self.assertIsInstance(chr1[2], str)  # GenBank accession
        self.assertIsInstance(chr1[3], str)  # RefSeq accession

    def test_hg38_chromosomes_structure(self) -> None:
        """Test hg38 chromosome data structure."""
        self.assertEqual(len(hg38_chromosomes), 25)  # 22 autosomes + X, Y, M

        chr1 = hg38_chromosomes[0]
        self.assertEqual(chr1[0], "1")  # chromosome name
        self.assertEqual(chr1[1], 248956422)  # length
        self.assertIsInstance(chr1[2], str)  # GenBank accession
        self.assertIsInstance(chr1[3], str)  # RefSeq accession

    def test_chromosome_names_consistent(self) -> None:
        """Test that chromosome names are consistent between hg19 and hg38."""
        hg19_names = [chrom[0] for chrom in hg19_chromosomes]
        hg38_names = [chrom[0] for chrom in hg38_chromosomes]
        self.assertEqual(hg19_names, hg38_names)


if __name__ == "__main__":
    unittest.main()
