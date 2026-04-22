# Integration tests for scripts/differential_wps.py.
#
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
import logging
import tempfile
import unittest

import numpy as np
import numpy.typing as npt
import pandas as pd

from pathlib import Path
from scripts.differential_wps import (
    load_wps_metric_matrix, perform_differential_analysis
)


def _write_sample_metrics(
    path: Path, regions: list[str], values: npt.NDArray[np.float64],
    extra_cols: dict[str, float] | None = None,
) -> None:
    """
    Write a minimal wps_metrics CSV (one row per region) with the requested
    metric values.
    """
    extras = extra_cols or {
        "n_peaks": 10, "median_peak_dist": 170.0,
        "median_prom": 0.1, "median_width": 40.0, "n_positions": 1000,
    }
    rows = []
    for region, v in zip(regions, values):
        rows.append({
            "region": region,
            "power_170": float(v),
            **extras,
        })
    pd.DataFrame(rows).to_csv(path, index=False)


class TestDifferentialWPS(unittest.TestCase):
    """End-to-end behaviour of differential_wps.py on synthetic metrics."""

    _logger: logging.Logger = logging.getLogger("test_differential_wps")

    def _build_cohort(
        self,
        out_dir: Path,
        n_group_a: int,
        n_group_b: int,
        diff_regions: list[str],
        null_regions: list[str],
        mean_a: float = 1.0,
        mean_b: float = 4.0,
        noise_sd: float = 0.2,
        seed: int = 0,
    ) -> tuple[list[str], list[str]]:
        """Write ``n_group_a`` + ``n_group_b`` per-sample CSVs to ``out_dir``.

        ``diff_regions`` carry mean ``mean_a`` in group A and ``mean_b``
        in group B; ``null_regions`` carry the same mean in both groups.
        """
        rng = np.random.default_rng(seed)
        all_regions = diff_regions + null_regions
        group_a_paths: list[str] = []
        group_b_paths: list[str] = []
        for i in range(n_group_a):
            vals = np.concatenate([
                rng.normal(mean_a, noise_sd, size=len(diff_regions)),
                rng.normal((mean_a + mean_b) / 2, noise_sd,
                           size=len(null_regions)),
            ])
            p = out_dir / f"wps_metrics_A{i}.csv"
            _write_sample_metrics(p, all_regions, vals)
            group_a_paths.append(str(p))
        for i in range(n_group_b):
            vals = np.concatenate([
                rng.normal(mean_b, noise_sd, size=len(diff_regions)),
                rng.normal((mean_a + mean_b) / 2, noise_sd,
                           size=len(null_regions)),
            ])
            p = out_dir / f"wps_metrics_B{i}.csv"
            _write_sample_metrics(p, all_regions, vals)
            group_b_paths.append(str(p))
        return group_a_paths, group_b_paths

    def test_load_wps_metric_matrix_shape_and_values(self) -> None:
        """load_wps_metric_matrix returns samples × regions with NaN where
        a region is missing from a sample."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            regions = ["r1", "r2", "r3"]
            _write_sample_metrics(
                tmp_p / "wps_metrics_A.csv", regions, np.array([1.0, 2.0, 3.0])
            )
            _write_sample_metrics(
                tmp_p / "wps_metrics_B.csv", regions[:2], np.array([4.0, 5.0])
            )
            mat = load_wps_metric_matrix(
                [str(tmp_p / "wps_metrics_A.csv"),
                 str(tmp_p / "wps_metrics_B.csv")],
                metric="power_170", group_name="test",
                logger=self._logger,
            )
            self.assertEqual(mat.shape, (2, 3))
            self.assertTrue(
                np.isnan(mat.loc["wps_metrics_B", "r3"])  # type: ignore
            )
            self.assertEqual(mat.loc["wps_metrics_A", "r1"], 1.0)
            self.assertEqual(mat.loc["wps_metrics_B", "r2"], 5.0)

    def test_differential_analysis_detects_known_signal(self) -> None:
        """perform_differential_analysis flags the seeded differential
        regions as significant and leaves null regions unflagged."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            diff_regions = [f"d{i}" for i in range(5)]
            null_regions = [f"n{i}" for i in range(20)]
            a_paths, b_paths = self._build_cohort(
                tmp_p, n_group_a=8, n_group_b=8,
                diff_regions=diff_regions, null_regions=null_regions,
            )

            mat_a = load_wps_metric_matrix(
                a_paths, "power_170", "A", self._logger
            )
            mat_b = load_wps_metric_matrix(
                b_paths, "power_170", "B", self._logger
            )
            results = perform_differential_analysis(
                mat_a, mat_b, metric="power_170",
                min_samples_per_group=3, logger=self._logger,
            )

            self.assertEqual(
                len(results), len(diff_regions) + len(null_regions)
            )
            self.assertIn("p_adjusted", results.columns)
            self.assertIn("effect_size", results.columns)
            self.assertIn("log_fold_change", results.columns)

            sig = set(results.loc[results["significant"], "region"])
            # All 5 differential regions should be significant.
            for r in diff_regions:
                self.assertIn(
                    r, sig,
                    f"expected seeded-differential region '{r}' to be sig"
                )
            # Null regions should overwhelmingly not be significant — allow
            # up to one false positive across 20 (FDR 5%).
            null_sig_fraction = sum(
                1 for r in null_regions if r in sig
            ) / len(null_regions)
            self.assertLessEqual(
                null_sig_fraction, 0.10,
                "too many null regions were flagged significant"
            )

            # log_fold_change should be positive for all differential
            # regions (group B > group A in our synthetic).
            d_rows = results.set_index("region").loc[diff_regions]
            self.assertTrue((d_rows["log_fold_change"] > 0).all())

    def test_min_samples_per_group_filter(self) -> None:
        """Regions with fewer than ``min_samples_per_group`` non-NaN values
        in either group are dropped."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            regions = ["full", "half_missing"]
            for i in range(5):
                _write_sample_metrics(
                    tmp_p / f"wps_metrics_A{i}.csv", regions,
                    np.array([1.0, 2.0]),
                )
                _write_sample_metrics(
                    tmp_p / f"wps_metrics_B{i}.csv", regions,
                    np.array([4.0, 5.0]),
                )
            # Re-write A4 so that ``half_missing`` is absent — four group-A
            # samples will now carry that region instead of five.
            _write_sample_metrics(
                tmp_p / "wps_metrics_A4.csv", regions[:1],
                np.array([1.0]),
            )
            a_paths = sorted(
                str(tmp_p / f"wps_metrics_A{i}.csv") for i in range(5)
            )
            b_paths = sorted(
                str(tmp_p / f"wps_metrics_B{i}.csv") for i in range(5)
            )

            mat_a = load_wps_metric_matrix(
                a_paths, "power_170", "A", self._logger
            )
            mat_b = load_wps_metric_matrix(
                b_paths, "power_170", "B", self._logger
            )

            # With min=4, both regions clear the filter (A has 5 and 4
            # non-NaN values, B has 5 and 5).
            results_lo = perform_differential_analysis(
                mat_a, mat_b, metric="power_170",
                min_samples_per_group=4, logger=self._logger,
            )
            self.assertEqual(
                sorted(results_lo["region"]), ["full", "half_missing"]
            )

            # Raising min to 5 drops ``half_missing`` (only 4 non-NaN in A)
            # while keeping ``full``.
            results_hi = perform_differential_analysis(
                mat_a, mat_b, metric="power_170",
                min_samples_per_group=5, logger=self._logger,
            )
            self.assertEqual(list(results_hi["region"]), ["full"])
