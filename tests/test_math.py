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
import json

import numpy as np
import numpy.typing as npt

from unittest.mock import patch, MagicMock
from pyfraglib.math import gaussian_mixture, mixture_cdf, \
                           mixture_cdf_wrapper, negative_log_likelihood, \
                           read_gmm_config, fit_gmm, \
                           jensen_shannon_divergence, goodness_of_fit_stats, \
                           plot_gmm, hessian, LARGE_DATASET_THRESHOLD


class TestGaussianMixture(unittest.TestCase):
    """Test Gaussian mixture model functions."""

    def test_gaussian_mixture_single_component(self) -> None:
        """Test Gaussian mixture with single component."""
        data: npt.NDArray[np.float64] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        params: list[float] = [3.0, 1.0, 1.0]
        n: int = 1
        pdf = gaussian_mixture(params, n, data)

        self.assertEqual(len(pdf), len(data))
        self.assertTrue(np.all(pdf > 0))
        peak_idx: np.int64 = np.argmax(pdf)
        self.assertEqual(data[peak_idx], 3.0)  # type: ignore

    def test_gaussian_mixture_multiple_components(self) -> None:
        """Test Gaussian mixture with multiple components."""
        data: npt.NDArray[np.float64] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        params: list[float] = [2.0, 4.0, 0.5, 0.5, 0.6, 0.4]
        n: int = 2

        pdf = gaussian_mixture(params, n, data)

        self.assertEqual(len(pdf), len(data))
        self.assertTrue(np.all(pdf > 0))

    def test_gaussian_mixture_assertions(self) -> None:
        """Test Gaussian mixture function assertions."""
        data: npt.NDArray[np.float64] = np.array([1.0, 2.0, 3.0])

        with self.assertRaises(AssertionError):
            gaussian_mixture([1.0, 1.0, 1.0], 0, data)

        with self.assertRaises(AssertionError):
            gaussian_mixture([1.0, 1.0], 1, data)

    def test_mixture_cdf(self) -> None:
        """Test mixture CDF calculation."""
        params: list[float] = [2.0, 1.0, 1.0]
        n: int = 1

        cdf_at_mean = mixture_cdf(2.0, params, n)
        self.assertAlmostEqual(cdf_at_mean, 0.5, places=2)
        x_values: list[float] = [0.0, 1.0, 2.0, 3.0, 4.0]
        cdf_values = [mixture_cdf(x, params, n) for x in x_values]

        for i in range(1, len(cdf_values)):
            self.assertGreaterEqual(cdf_values[i], cdf_values[i-1])

    def test_mixture_cdf_wrapper(self) -> None:
        """Test mixture CDF wrapper function."""
        params: list[float] = [2.0, 1.0, 1.0]
        n: int = 1
        data: list[float] = [1.0, 2.0, 3.0]
        cdf_array: npt.NDArray[np.float64] = \
            mixture_cdf_wrapper(data, params, n)

        self.assertEqual(len(cdf_array), len(data))

        for i, x in enumerate(data):
            expected: float = mixture_cdf(x, params, n)
            self.assertAlmostEqual(
                cdf_array[i], expected, places=6  # type: ignore
            )

    def test_negative_log_likelihood(self) -> None:
        """Test negative log likelihood calculation."""
        np.random.seed(42)
        data: npt.NDArray[np.float64] = np.random.normal(2.0, 1.0, 100)
        params: list[float] = [2.0, 1.0, 1.0]
        n: int = 1
        norm_const: float = 1.0

        nll = negative_log_likelihood(params, n, data, norm_const)
        self.assertIsInstance(nll, float)
        self.assertGreater(nll, 0)

        bad_params: list[float] = [0.0, 0.1, 1.0]
        bad_nll = negative_log_likelihood(bad_params, n, data, norm_const)
        self.assertGreater(bad_nll, nll)


class TestGMMConfig(unittest.TestCase):
    """Test GMM configuration reading."""

    def create_test_config(self, config_data: dict[str, object]) -> str:
        """Create a temporary config file for testing."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            json.dump(config_data, f)
            return f.name

    def test_read_gmm_config_valid(self) -> None:
        """Test reading valid GMM configuration."""
        config_data = {
            "number_of_gaussians": 2,
            "subsample_percentage": 0.1,
            "means_lower_bounds": [100.0, 200.0],
            "means_upper_bounds": [150.0, 250.0],
            "std_lower_bounds": [10.0, 20.0],
            "std_upper_bounds": [30.0, 40.0],
            "initial_means": [120.0, 220.0],
            "initial_pis": [0.6, 0.4]
        }

        config_path = self.create_test_config(config_data)
        try:
            n, ssp, initials, bounds = read_gmm_config(config_path)

            self.assertEqual(n, 2)
            self.assertEqual(ssp, 0.1)
            self.assertEqual(len(initials), 6)  # 2*3 parameters
            self.assertEqual(len(bounds), 6)  # 2*3 bounds
            self.assertEqual(initials[0], 120.0)  # First mean
            self.assertEqual(initials[1], 220.0)  # Second mean
            self.assertEqual(initials[4], 0.6)    # First pi
            self.assertEqual(initials[5], 0.4)    # Second pi
            self.assertEqual(bounds[0], (100.0, 150.0))  # First mean bounds
            self.assertEqual(bounds[1], (200.0, 250.0))  # Second mean bounds
            self.assertEqual(bounds[4], (0.0, 1.0))      # Pi bounds

        finally:
            os.unlink(config_path)

    @patch("pyfraglib.math.fail")
    def test_read_gmm_config_invalid_gaussians(
        self, mock_fail: MagicMock
    ) -> None:
        """Test invalid number of Gaussians."""
        mock_fail.side_effect = SystemExit(1)
        config_data: dict[str, object] = {"number_of_gaussians": 0}
        config_path = self.create_test_config(config_data)
        try:
            with self.assertRaises(SystemExit):
                read_gmm_config(config_path)
        finally:
            os.unlink(config_path)

    @patch("pyfraglib.math.fail")
    def test_read_gmm_config_invalid_subsample(
        self, mock_fail: MagicMock
    ) -> None:
        """Test invalid subsample percentage."""
        mock_fail.side_effect = SystemExit(1)
        config_data: dict[str, object] = {
            "number_of_gaussians": 1,
            "subsample_percentage": 1.5
        }
        config_path = self.create_test_config(config_data)
        try:
            with self.assertRaises(SystemExit):
                read_gmm_config(config_path)
        finally:
            os.unlink(config_path)

    @patch("pyfraglib.math.fail")
    def test_read_gmm_config_missing_fields(
        self, mock_fail: MagicMock
    ) -> None:
        """Test missing required fields."""
        mock_fail.side_effect = SystemExit(1)
        config_data: dict[str, object] = {
            "number_of_gaussians": 1,
            "subsample_percentage": 0.1
        }
        config_path = self.create_test_config(config_data)
        try:
            with self.assertRaises(SystemExit):
                read_gmm_config(config_path)
        finally:
            os.unlink(config_path)

    @patch("pyfraglib.math.fail")
    def test_read_gmm_config_wrong_array_length(
        self, mock_fail: MagicMock
    ) -> None:
        """Test arrays with wrong length."""
        mock_fail.side_effect = SystemExit(1)
        config_data: dict[str, object] = {
            "number_of_gaussians": 2,
            "subsample_percentage": 0.1,
            "means_lower_bounds": [100.0],
            "means_upper_bounds": [150.0, 250.0],
            "std_lower_bounds": [10.0, 20.0],
            "std_upper_bounds": [30.0, 40.0],
            "initial_means": [120.0, 220.0],
            "initial_pis": [0.6, 0.4]
        }
        config_path = self.create_test_config(config_data)
        try:
            with self.assertRaises(SystemExit):
                read_gmm_config(config_path)
        finally:
            os.unlink(config_path)


class TestGMMFitting(unittest.TestCase):
    """Test GMM fitting functionality."""

    def create_test_config_file(self) -> str:
        """Create a valid test configuration file."""
        config_data = {
            "number_of_gaussians": 1,
            "subsample_percentage": 1.0,  # Use all data for testing
            "means_lower_bounds": [140.0],
            "means_upper_bounds": [180.0],
            "std_lower_bounds": [5.0],
            "std_upper_bounds": [25.0],
            "initial_means": [160.0],
            "initial_pis": [1.0]
        }

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as f:
            json.dump(config_data, f)
            return f.name

    @patch("pyfraglib.math.minimize")
    def test_fit_gmm_basic(self, mock_minimize: MagicMock) -> None:
        """Test basic GMM fitting functionality."""
        np.random.seed(42)

        mock_result = MagicMock()
        mock_result.x = [160.0, 15.0, 1.0]  # type: ignore
        mock_minimize.return_value = mock_result
        data: npt.NDArray[np.float64] = np.random.normal(160.0, 15.0, 100)
        config_path = self.create_test_config_file()
        try:
            result, n, params, data_sample = fit_gmm(data, config_path)
            self.assertEqual(n, 1)
            self.assertEqual(len(params), 3)  # mean, std, pi
            self.assertEqual(len(data_sample), 100)  # all data is used
            mock_minimize.assert_called_once()
        finally:
            os.unlink(config_path)

    def test_hessian(self) -> None:
        """Test hessian function (returns zeros)."""
        params = [1.0, 2.0, 3.0]
        n = 1
        data = [1.0, 2.0, 3.0]
        hess = hessian(params, n, data)
        self.assertEqual(hess.shape, (3, 3))
        self.assertTrue(np.all(hess == 0))  # type: ignore


class TestDivergenceMetrics(unittest.TestCase):
    """Test divergence and distance metrics."""

    def test_jensen_shannon_divergence_identical(self) -> None:
        """Test JS divergence between identical distributions."""
        p: npt.NDArray[np.float64] = np.array([0.5, 0.3, 0.2])
        q: npt.NDArray[np.float64] = np.array([0.5, 0.3, 0.2])
        js_div = jensen_shannon_divergence(p, q)
        self.assertAlmostEqual(js_div, 0.0, places=10)

    def test_jensen_shannon_divergence_different(self) -> None:
        """Test JS divergence between different distributions."""
        p: npt.NDArray[np.float64] = np.array([0.8, 0.1, 0.1])
        q: npt.NDArray[np.float64] = np.array([0.1, 0.8, 0.1])
        js_div = jensen_shannon_divergence(p, q)
        self.assertGreater(js_div, 0.0)
        self.assertLessEqual(js_div, np.log(2))  # type: ignore

    def test_jensen_shannon_divergence_normalization(self) -> None:
        """Test JS divergence with unnormalized distributions."""
        p: npt.NDArray[np.float64] = np.array([2.0, 1.0, 1.0])
        q: npt.NDArray[np.float64] = np.array([1.0, 2.0, 1.0])
        js_div = jensen_shannon_divergence(p, q)
        self.assertGreater(js_div, 0.0)
        self.assertLessEqual(js_div, np.log(2))  # type: ignore


class TestGoodnessOfFit(unittest.TestCase):
    """Test goodness of fit statistics."""

    def test_goodness_of_fit_stats(self) -> None:
        """Test goodness of fit statistics calculation."""
        np.random.seed(42)
        data: npt.NDArray[np.float64] = np.random.normal(150.0, 20.0, 1000)
        params: list[float] = [150.0, 20.0, 1.0]
        n: int = 1
        stats: dict[str, object] = goodness_of_fit_stats(data, params, n)
        expected_keys: list[str] = [
            "kolmogorov_smirnov_statistic",
            "kolmogorov_smirnov_p_value",
            "wasserstein_distance",
            "jensen_shannon_divergence"
        ]

        for key in expected_keys:
            self.assertIn(key, stats)
            self.assertIsInstance(stats[key], (int, float))

        self.assertGreater(
            stats["kolmogorov_smirnov_p_value"], 0.01  # type: ignore
        )
        self.assertLess(stats["wasserstein_distance"], 10.0)  # type: ignore
        self.assertLess(
            stats["jensen_shannon_divergence"], 0.15  # type: ignore
        )


class TestPlotting(unittest.TestCase):
    """Test plotting functions."""

    @patch("matplotlib.pyplot.savefig")  # type: ignore
    @patch("matplotlib.pyplot.figure")
    def test_plot_gmm(
        self, mock_figure: MagicMock, mock_savefig: MagicMock
    ) -> None:
        """Test GMM plotting function."""
        np.random.seed(42)

        mock_fig = MagicMock()
        mock_fig.dpi = 100
        mock_figure.return_value = mock_fig
        data: npt.NDArray[np.float64] = np.random.normal(150.0, 20.0, 100)
        params: list[float] = [150.0, 20.0, 1.0]
        n: int = 1

        with tempfile.TemporaryDirectory() as temp_dir:
            plot_gmm(data, n, params, temp_dir, "test_sample")
            expected_path: str = os.path.join(
                temp_dir, "test_sample_gmm_frags_len.png"
            )
            mock_fig.savefig.assert_called_once_with(  # type: ignore
                expected_path, dpi=100
            )

    @patch("matplotlib.pyplot.savefig")  # type: ignore
    @patch("matplotlib.pyplot.figure")
    def test_plot_gmm_multiple_components(
        self, mock_figure: MagicMock, mock_savefig: MagicMock
    ) -> None:
        """Test GMM plotting with multiple components."""
        np.random.seed(42)

        mock_fig = MagicMock()
        mock_fig.dpi = 100
        mock_figure.return_value = mock_fig

        data1: npt.NDArray[np.float64] = np.random.normal(120.0, 15.0, 50)
        data2: npt.NDArray[np.float64] = np.random.normal(180.0, 20.0, 50)
        data: npt.NDArray[np.float64] = np.concatenate(
            [data1, data2]  # type: ignore
        )
        params: list[float] = [120.0, 180.0, 15.0, 20.0, 0.5, 0.5]
        n: int = 2

        with tempfile.TemporaryDirectory() as temp_dir:
            plot_gmm(data, n, params, temp_dir, "bimodal_test")
            mock_fig.savefig.assert_called_once()  # type: ignore


class TestConstants(unittest.TestCase):
    """Test module constants."""

    def test_large_dataset_threshold(self) -> None:
        """Test large dataset threshold constant."""
        self.assertEqual(LARGE_DATASET_THRESHOLD, 1_000_000)
        self.assertIsInstance(LARGE_DATASET_THRESHOLD, int)


if __name__ == "__main__":
    unittest.main()
