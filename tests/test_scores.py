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
import pysam

import pandas as pd
import numpy as np

from unittest.mock import Mock, patch, MagicMock
from tests.test_fixtures import create_mock_fragment
from pyfraglib.fragment import FragmentList
from pyfraglib.scores import motif_diversity, windowed_protection_score, \
                             windowed_protection_score_fast


class TestMotifDiversity(unittest.TestCase):
    """Test motif diversity scoring functions."""

    def create_fragment_list_with_motifs(self) -> FragmentList:
        """Create a FragmentList with specific end motifs for testing."""
        fragment_list = FragmentList()
        motifs_5p = [
            "ATCGAAAA", "GCTTTTTT", "TAGACCCC", "CGATGGGG", "AAAATTTT"
        ]
        motifs_3p = [
            "AAAACGAT", "TTTTGCTA", "CCCCTAGA", "GGGGCGAT", "TTTTAAAA"
        ]

        for i, (motif5p, motif3p) in enumerate(zip(motifs_5p, motifs_3p)):
            fragment = create_mock_fragment(
                start_pos=100 + i * 100,
                end_pos=200 + i * 100,
                chrom="1",
                end5p=motif5p,
                end3p=motif3p,
                is_bogus=False
            )
            fragment_list.append(fragment)

        return fragment_list

    def test_motif_diversity_shannon(self) -> None:
        """Test motif diversity calculation with Shannon entropy."""
        fragment_list = self.create_fragment_list_with_motifs()
        index_5p, index_3p = motif_diversity(
            fragment_list, "test_sample", "shannon"
        )

        self.assertIsInstance(index_5p, float)
        self.assertIsInstance(index_3p, float)
        self.assertGreater(index_5p, 0.0)
        self.assertGreater(index_3p, 0.0)
        self.assertLess(index_5p, 2.0)
        self.assertLess(index_3p, 2.0)

    def test_motif_diversity_simpson(self) -> None:
        """Test motif diversity calculation with Simpson index."""
        fragment_list = self.create_fragment_list_with_motifs()
        index_5p, index_3p = motif_diversity(
            fragment_list, "test_sample", "simpson"
        )

        self.assertIsInstance(index_5p, float)
        self.assertIsInstance(index_3p, float)
        self.assertGreater(index_5p, 0.0)
        self.assertGreater(index_3p, 0.0)
        self.assertLessEqual(index_5p, 1.0)
        self.assertLessEqual(index_3p, 1.0)

    def test_motif_diversity_uniform_distribution(self) -> None:
        """Test motif diversity with uniform motif distribution."""
        fragment_list = FragmentList()
        motifs = ["AAAA", "TTTT", "CCCC", "GGGG"]
        for i in range(16):
            motif = motifs[i % 4]
            fragment = create_mock_fragment(
                start_pos=100 + i * 100,
                end_pos=200 + i * 100,
                end5p=motif,
                end3p=motif,
                is_bogus=False
            )
            fragment_list.append(fragment)

        shannon_5p, shannon_3p = motif_diversity(
            fragment_list, "uniform", "shannon"
        )
        simpson_5p, simpson_3p = motif_diversity(
            fragment_list, "uniform", "simpson"
        )
        self.assertGreater(shannon_5p, 1.0)
        self.assertLess(simpson_5p, 0.5)

    def test_motif_diversity_single_motif(self) -> None:
        """Test motif diversity with single dominant motif."""
        fragment_list = FragmentList()
        for i in range(10):
            fragment = create_mock_fragment(
                start_pos=100 + i * 100,
                end_pos=200 + i * 100,
                end5p="ATCGTTTT",
                end3p="TTTTCGTA",
                is_bogus=False
            )
            fragment_list.append(fragment)

        shannon_5p, shannon_3p = motif_diversity(
            fragment_list, "single", "shannon"
        )
        simpson_5p, simpson_3p = motif_diversity(
            fragment_list, "single", "simpson"
        )
        self.assertAlmostEqual(shannon_5p, 0.0, places=10)
        self.assertAlmostEqual(simpson_5p, 1.0, places=10)

    @patch("pyfraglib.scores.fail")
    def test_motif_diversity_invalid_index(self, mock_fail: MagicMock) -> None:
        """Test motif diversity with invalid index function."""
        mock_fail.side_effect = SystemExit(1)
        fragment_list = self.create_fragment_list_with_motifs()

        with self.assertRaises(SystemExit):
            motif_diversity(fragment_list, "test", "invalid_index")

    def test_motif_diversity_with_bogus_fragments(self) -> None:
        """Test bogus fragment exclusion from diversity calculation."""
        fragment_list = FragmentList()
        normal_fragment = create_mock_fragment(
            end5p="ATCGTTTT", end3p="TTTTCGTA", is_bogus=False
        )
        fragment_list.append(normal_fragment)

        bogus_fragment = create_mock_fragment(
            end5p="NNNNAAAA", end3p="AAAANNNN", is_bogus=True
        )
        fragment_list.append(bogus_fragment)

        shannon_5p, shannon_3p = motif_diversity(
            fragment_list, "test", "shannon"
        )

        self.assertAlmostEqual(shannon_5p, 0.0, places=10)
        self.assertAlmostEqual(shannon_3p, 0.0, places=10)


class TestWindowedProtectionScore(unittest.TestCase):
    """Test windowed protection score functions."""

    def create_mock_tabix_file(self) -> Mock:
        """Create a mock TabixFile for testing."""
        mock_tabix: Mock = Mock(spec=pysam.TabixFile)  # type: ignore
        regions = [
            "1\t1000\t2000\tregion1",
            "1\t3000\t4000\tregion2",
            "2\t1000\t2000\tregion3"
        ]
        mock_tabix.fetch.return_value = regions  # type: ignore
        return mock_tabix

    def create_test_fragments_for_wps(self) -> FragmentList:
        """Create fragments specifically for WPS testing."""
        fragment_list = FragmentList()
        positions = [
            (950, 1050),
            (1500, 1700),
            (900, 2100),
            (1800, 1950),
        ]

        for i, (start, end) in enumerate(positions):
            fragment = create_mock_fragment(
                start_pos=start,
                end_pos=end,
                chrom="1",
                length=end - start,
                is_bogus=False
            )
            fragment_list.append(fragment)

        return fragment_list

    @patch("pyfraglib.scores.chromosome_maps_to_df")  # type: ignore
    @patch("pyfraglib.scores.create_chromosome_map")
    def test_windowed_protection_score_fast(
        self, mock_create_map: Mock, mock_maps_to_df: Mock
    ) -> None:
        """Test fast windowed protection score calculation."""
        mock_map = np.zeros(250000000, dtype=np.int64)
        mock_create_map.return_value = \
            {"1": mock_map, "2": mock_map}  # type: ignore
        mock_df: pd.DataFrame = pd.DataFrame({
            "chrom": ["1", "1"],
            "pos": [1500, 3500],
            "wps": [1.5, -0.5]
        })
        mock_maps_to_df.return_value = mock_df

        fragment_list = self.create_test_fragments_for_wps()
        mock_regions = self.create_mock_tabix_file()
        _ = windowed_protection_score_fast(
            fragment_list, mock_regions, win_size=120, genome="hg19"
        )
        self.assertEqual(mock_create_map.call_count, 4)
        mock_maps_to_df.assert_called_once()

    def test_windowed_protection_score_dispatch(self) -> None:
        """Test dispatch to fast implementation."""
        fragment_list = self.create_test_fragments_for_wps()
        mock_regions = self.create_mock_tabix_file()

        with patch(
            "pyfraglib.scores.windowed_protection_score_fast"
        ) as mock_fast:
            mock_fast.return_value = pd.DataFrame()
            _ = windowed_protection_score(fragment_list, mock_regions)
            mock_fast.assert_called_once_with(
                fragment_list, mock_regions, 120, "hg19"
            )

    def test_windowed_protection_score_custom_params(self) -> None:
        """Test windowed protection score with custom parameters."""
        fragment_list = self.create_test_fragments_for_wps()
        mock_regions = self.create_mock_tabix_file()
        with patch(
            "pyfraglib.scores.windowed_protection_score_fast"
        ) as mock_fast:
            mock_fast.return_value = pd.DataFrame()
            windowed_protection_score(
                fragment_list, mock_regions, win_size=240, genome="hg38"
            )
            mock_fast.assert_called_once_with(
                fragment_list, mock_regions, 240, "hg38"
            )

    @patch("pyfraglib.scores.get_chromosome_length")
    def test_windowed_protection_score_fragment_processing(
        self, mock_get_length: Mock
    ) -> None:
        """Test fragment processing logic in WPS calculation."""
        mock_get_length.return_value = 250000000
        fragment_list = FragmentList()
        short_fragment = create_mock_fragment(
            start_pos=1000, end_pos=1100, chrom="1", length=100, is_bogus=False
        )
        long_fragment = create_mock_fragment(
            start_pos=2000, end_pos=2200, chrom="1", length=200, is_bogus=False
        )
        bogus_fragment = create_mock_fragment(
            start_pos=3000, end_pos=3100, chrom="1", length=100, is_bogus=True
        )
        fragment_list.append(short_fragment)
        fragment_list.append(long_fragment)
        fragment_list.append(bogus_fragment)

        with patch(
            "pyfraglib.scores.create_chromosome_map"
        ) as mock_create_map:
            mock_map = np.zeros(250000000, dtype=np.int64)
            mock_create_map.return_value = {"1": mock_map}  # type: ignore

            with patch(
                "pyfraglib.scores.chromosome_maps_to_df"
            ) as mock_maps_to_df:
                mock_maps_to_df.return_value = pd.DataFrame()
                mock_regions = self.create_mock_tabix_file()
                windowed_protection_score_fast(fragment_list, mock_regions)
                self.assertEqual(mock_create_map.call_count, 4)


class TestWPSHelperFunctions(unittest.TestCase):
    """Test helper functions for WPS calculation."""

    def test_wps_assertion_win_size(self) -> None:
        """Test that WPS functions assert valid window size."""
        fragment_list = FragmentList()
        mock_regions = Mock()

        with self.assertRaises(AssertionError):
            windowed_protection_score_fast(
                fragment_list, mock_regions, win_size=0
            )

        with self.assertRaises(AssertionError):
            windowed_protection_score_fast(
                fragment_list, mock_regions, win_size=-10
            )


class TestWPSIntegration(unittest.TestCase):
    """Test WPS integration scenarios."""

    def create_test_fragments_for_wps(self) -> FragmentList:
        """Create fragments specifically for WPS testing."""
        fragment_list: FragmentList = FragmentList()
        positions = [
            (950, 1050),
            (1500, 1700),
            (900, 2100),
            (1800, 1950),
        ]

        for i, (start, end) in enumerate(positions):
            fragment = create_mock_fragment(
                start_pos=start,
                end_pos=end,
                chrom="1",
                length=end - start,
                is_bogus=False
            )
            fragment_list.append(fragment)

        return fragment_list

    def test_wps_empty_fragment_list(self) -> None:
        """Test WPS calculation with empty fragment list."""
        fragment_list = FragmentList()
        mock_regions = Mock(spec=pysam.TabixFile)  # type: ignore
        mock_regions.fetch.return_value = (  # type: ignore
            ["1\t1000\t2000\ttest"]  # type: ignore
        )

        with patch(
            "pyfraglib.scores.create_chromosome_map"
        ) as mock_create_map:
            mock_map = np.zeros(250000000, dtype=np.int64)
            mock_create_map.return_value = {"1": mock_map}  # type: ignore

            with patch(
                "pyfraglib.scores.chromosome_maps_to_df"
            ) as mock_maps_to_df:
                mock_maps_to_df.return_value = pd.DataFrame()
                _ = windowed_protection_score_fast(fragment_list, mock_regions)

    def test_wps_different_genomes(self) -> None:
        """Test WPS calculation with different genome versions."""
        fragment_list = self.create_test_fragments_for_wps()
        mock_regions = Mock(spec=pysam.TabixFile)  # type: ignore
        mock_regions.fetch.return_value = (  # type: ignore
            ["1\t1000\t2000\ttest"]  # type: ignore
        )

        for genome in ["hg19", "hg38"]:
            with patch(
                "pyfraglib.scores.create_chromosome_map"
            ) as mock_create_map:
                mock_map = np.zeros(250000000, dtype=np.int64)
                mock_create_map.return_value = {"1": mock_map}  # type: ignore

                with patch(
                    "pyfraglib.scores.chromosome_maps_to_df"
                ) as mock_maps_to_df:
                    mock_maps_to_df.return_value = pd.DataFrame()
                    _ = windowed_protection_score_fast(
                        fragment_list, mock_regions, genome=genome
                    )


if __name__ == "__main__":
    unittest.main()
