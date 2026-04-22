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
import unittest

import numpy as np
import numpy.typing as npt
import pandas as pd

from sklearn.decomposition import NMF  # type: ignore

from scripts.nmf_fragment_lengths import (
    perform_nmf_analysis, project_onto_signatures,
)


class TestNMFMixtureFractionRecovery(unittest.TestCase):
    """
    Integration tests for perform_nmf_analysis mixture fraction recovery.

    This test uses two perfectly non-overlapping synthetic fragment length
    distributions (windows at 100-119 bp and 200-219 bp).  With no shared
    features, NMF must recover the ground-truth mixture fractions exactly.
    This demonstrates that the implementation is correct and that the
    ~10 pp overestimation observed on real cfDNA data is caused by the
    substantial overlap between the hematopoietic and tumor-like fragment
    length distributions.
    """

    _logger: logging.Logger = logging.getLogger("test_nmf")

    def _make_data(
        self, n_pure_a: int, n_pure_b: int, n_mixed: int, frac_b: float,
        counts: int = 10_000
    ) -> pd.DataFrame:
        """
        Build a (n_pure_a + n_pure_b + n_mixed) x 120 fragment-length count
        matrix.

        Signature A: uniform counts over lengths 100-119 bp.
        Signature B: uniform counts over lengths 200-219 bp.
        The two windows are completely non-overlapping.

        Pure-A samples contain only Sig A.
        Pure-B samples contain only Sig B.
        Mixed samples contain (1 - frac_b) Sig A and frac_b Sig B.

        Both pure-A and pure-B anchor samples are required so that the NMF
        separability condition is satisfied and the decomposition is unique.
        """
        lengths = list(range(100, 220))
        n_samples = n_pure_a + n_pure_b + n_mixed
        data = np.zeros((n_samples, len(lengths)), dtype=np.float64)

        sig_a_idx = [i for i, ln in enumerate(lengths) if 100 <= ln < 120]
        sig_b_idx = [i for i, ln in enumerate(lengths) if 200 <= ln < 220]

        for i in range(n_pure_a):
            data[i, sig_a_idx] = counts / len(sig_a_idx)

        for i in range(n_pure_a, n_pure_a + n_pure_b):
            data[i, sig_b_idx] = counts / len(sig_b_idx)

        for i in range(n_pure_a + n_pure_b, n_samples):
            data[i, sig_a_idx] = counts * (1.0 - frac_b) / len(sig_a_idx)
            data[i, sig_b_idx] = counts * frac_b / len(sig_b_idx)

        return pd.DataFrame(data, columns=[str(ln) for ln in lengths])

    def _h_weighted_fractions(
        self, W: npt.NDArray[np.float64], H: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """
        Return H-weighted row-normalised W (true mixture fractions).

        Each column of W is scaled by the corresponding H row sum before
        row-normalisation, so that differences in signature magnitude are
        accounted for.
        """
        H_sums: npt.NDArray[np.float64] = H.sum(axis=1)
        weighted: npt.NDArray[np.float64] = W * H_sums[np.newaxis, :]
        row_sums: npt.NDArray[np.float64] = weighted.sum(axis=1, keepdims=True)
        return weighted / np.where(row_sums > 0, row_sums, 1.0)

    def test_non_overlapping_signatures_recover_exact_fractions(self) -> None:
        """
        NMF must recover ground-truth mixture fractions to within 1 pp when
        the two component signatures share no features (non-overlapping
        windows). Both pure-A and pure-B anchor samples are included so that
        the NMF separability condition is satisfied and the decomposition is
        unique.
        """
        TRUE_FRAC_B = 0.15
        N_PURE_A = 10
        N_PURE_B = 10
        N_MIXED = 10
        TOLERANCE = 0.01

        data_df = self._make_data(N_PURE_A, N_PURE_B, N_MIXED, TRUE_FRAC_B)
        W, H, error = perform_nmf_analysis(
            data_df, n_components=2, logger=self._logger
        )
        W_frac = self._h_weighted_fractions(W, H)

        lengths = [int(c) for c in data_df.columns]
        sig_b_mask: npt.NDArray[np.bool_] = np.array(
            [200 <= ln < 220 for ln in lengths]
        )
        sig_b_comp = int(np.argmax(H[:, sig_b_mask].sum(axis=1)))

        pure_a_frac_b: npt.NDArray[np.float64] = W_frac[:N_PURE_A, sig_b_comp]
        pure_b_frac_b: npt.NDArray[np.float64] = W_frac[
            N_PURE_A:N_PURE_A + N_PURE_B, sig_b_comp
        ]
        mixed_frac_b: npt.NDArray[np.float64] = \
            W_frac[N_PURE_A + N_PURE_B:, sig_b_comp]

        self.assertTrue(
            np.all(pure_a_frac_b < TOLERANCE),
            f"A samples should have ~0% Sig B; got {pure_a_frac_b.round(4)}"
        )
        self.assertTrue(
            np.all(np.abs(pure_b_frac_b - 1.0) < TOLERANCE),
            f"B samples should have ~100% Sig B; got {pure_b_frac_b.round(4)}"
        )
        self.assertTrue(
            np.all(np.abs(mixed_frac_b - TRUE_FRAC_B) < TOLERANCE),
            f"Mixed samples should have ~{TRUE_FRAC_B*100:.0f}% Sig B; "
            f"got {mixed_frac_b.round(4)}"
        )
        self.assertLess(
            error, 1e-6,
            f"Reconstruction error should be ~0 for non-overlapping data; "
            f"got {error:.2e}"
        )

    def _make_overlapping_data(
        self,
        n_pure_a: int,
        n_mixed: int,
        frac_b: float,
        counts: int = 10_000,
    ) -> pd.DataFrame:
        """
        Build a (n_pure_a + n_mixed) x 301 fragment-length count matrix using
        Gaussian-shaped signatures that substantially overlap, mimicking the
        real cfDNA scenario (hematopoietic vs. tumor-like tissue profiles).

        Signature A (hematopoietic-like): Gaussian, mean=167 bp, std=10 bp.
        Signature B (tumor-like):         Gaussian, mean=150 bp, std=18 bp.

        Pure-A samples contain only Sig A; no pure-B samples are included,
        mirroring the real cohort design (no pure-tumor controls).
        """
        lengths = np.arange(100, 401)
        sig_a: npt.NDArray[np.float64] = np.exp(
            -0.5 * ((lengths - 167.0) / 10.0) ** 2
        )
        sig_b: npt.NDArray[np.float64] = np.exp(
            -0.5 * ((lengths - 150.0) / 18.0) ** 2
        )
        sig_a = sig_a / sig_a.sum() * counts
        sig_b = sig_b / sig_b.sum() * counts

        n_samples = n_pure_a + n_mixed
        data = np.zeros((n_samples, len(lengths)), dtype=np.float64)
        for i in range(n_pure_a):
            data[i] = sig_a
        for i in range(n_pure_a, n_samples):
            data[i] = (1.0 - frac_b) * sig_a + frac_b * sig_b

        return pd.DataFrame(data, columns=[str(ln) for ln in lengths])

    def test_overlapping_signatures_overestimate_fractions(self) -> None:
        """
        NMF overestimates minority-component fractions when the two component
        signatures overlap substantially and no pure minority-component anchor
        samples are present (i.e. the separability condition is violated). This
        mirrors real cohort designs (no pure tumor samples).
        """
        TRUE_FRAC_B = 0.15
        N_PURE_A = 10
        N_MIXED = 10
        MIN_OVERESTIMATION = 0.05
        MAX_FRACTION = 0.60

        data_df = self._make_overlapping_data(N_PURE_A, N_MIXED, TRUE_FRAC_B)
        W, H, _ = perform_nmf_analysis(
            data_df, n_components=2, logger=self._logger
        )
        W_frac = self._h_weighted_fractions(W, H)

        lengths_arr: npt.NDArray[np.float64] = np.array(
            [float(c) for c in data_df.columns]
        )
        h_peaks: npt.NDArray[np.float64] = np.array(
            [lengths_arr[int(np.argmax(H[k, :]))] for k in range(H.shape[0])]
        )
        sig_b_comp = int(np.argmin(np.abs(h_peaks - 150.0)))

        mixed_frac_b: npt.NDArray[np.float64] = W_frac[N_PURE_A:, sig_b_comp]
        mean_recovered = float(mixed_frac_b.mean())

        self.assertGreater(
            mean_recovered,
            TRUE_FRAC_B + MIN_OVERESTIMATION,
            f"Overlapping signatures without pure-B anchors should cause "
            f"overestimation (> {TRUE_FRAC_B + MIN_OVERESTIMATION:.2f}); "
            f"got mean={mean_recovered:.4f}"
        )
        self.assertLess(
            mean_recovered,
            MAX_FRACTION,
            f"Overestimation should be bounded (< {MAX_FRACTION:.2f}); "
            f"got mean={mean_recovered:.4f}"
        )


class TestNMFProjection(unittest.TestCase):
    """
    Equivalence tests for ``project_onto_signatures``.
    """

    _logger: logging.Logger = logging.getLogger("test_nmf")

    @staticmethod
    def _row_normalise(
        X: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        row_sums = X.sum(axis=1, keepdims=True)
        return X / np.where(row_sums > 0, row_sums, 1.0)

    def _fit_and_project_both(
        self,
        X_fit: npt.NDArray[np.float64],
        X_new: npt.NDArray[np.float64],
        n_components: int,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64],
               npt.NDArray[np.float64]]:
        X_fit_n = self._row_normalise(X_fit)
        X_new_n = self._row_normalise(X_new)
        model = NMF(
            n_components=n_components, init="nndsvd",
            random_state=0, max_iter=5000, tol=1e-6,
        )
        model.fit(X_fit_n)
        H = model.components_
        W_sklearn = model.transform(X_new_n)
        W_ours, _ = project_onto_signatures(
            pd.DataFrame(X_new), H, self._logger,
        )
        return W_ours, W_sklearn, H

    def test_projection_matches_sklearn_nonoverlapping(self) -> None:
        """Projection onto non-overlapping signatures recovers the
        ground-truth mixing fractions to within a tight tolerance, and is
        numerically equivalent to sklearn's ``NMF.transform``."""
        rng = np.random.default_rng(42)
        n_features = 40
        sig_a = np.zeros(n_features)
        sig_a[:20] = 1.0
        sig_b = np.zeros(n_features)
        sig_b[20:] = 1.0

        fit_rows = []
        for _ in range(5):
            fit_rows.append(sig_a * 1000 + rng.poisson(5, n_features))
            fit_rows.append(sig_b * 1000 + rng.poisson(5, n_features))
        X_fit = np.vstack(fit_rows).astype(float)

        true_frac_b = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        X_new = np.stack([
            (1 - f) * sig_a * 1000 + f * sig_b * 1000 +
            rng.poisson(5, n_features)
            for f in true_frac_b
        ]).astype(float)

        W_ours, W_sklearn, H = self._fit_and_project_both(
            X_fit, X_new, n_components=2,
        )

        # Normalise to mixture fractions for both.
        def _to_frac(W: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            scaled = W * H.sum(axis=1)[None, :]
            return scaled / scaled.sum(axis=1, keepdims=True)  # type: ignore

        f_ours = _to_frac(W_ours)
        f_sklearn = _to_frac(W_sklearn)

        # Align components by matching each to its nearest sklearn counterpart.
        sig_b_idx = int(np.argmin(np.abs(H @ sig_b - (H @ sig_b).max())))
        rec_b_ours = f_ours[:, sig_b_idx]
        rec_b_skl = f_sklearn[:, sig_b_idx]

        np.testing.assert_allclose(
            rec_b_ours, true_frac_b, atol=0.02,
            err_msg="scipy-NNLS projection should recover ground-truth "
                    "fractions to within 2 pp on non-overlapping signatures")
        np.testing.assert_allclose(
            rec_b_ours, rec_b_skl, atol=1e-4,
            err_msg="scipy-NNLS projection must match sklearn's NMF.transform")

    def test_projection_matches_sklearn_random_signatures(self) -> None:
        """On random dense signatures, the scipy-NNLS projection matches
        sklearn's ``NMF.transform`` to within numerical tolerance, and
        reconstructs the fit set with a small squared error."""
        rng = np.random.default_rng(7)
        n_components = 3
        n_features = 60
        H_true = rng.uniform(0.1, 1.0, size=(n_components, n_features))
        frac_fit = rng.dirichlet(np.ones(n_components), size=12)
        X_fit = \
            frac_fit @ H_true * 1000 + rng.poisson(3, size=(12, n_features))
        frac_new = rng.dirichlet(np.ones(n_components), size=6)
        X_new = frac_new @ H_true * 1000 + rng.poisson(3, size=(6, n_features))

        W_ours, W_sklearn, _H = self._fit_and_project_both(
            X_fit, X_new, n_components=n_components,
        )

        np.testing.assert_allclose(
            W_ours, W_sklearn, atol=5e-3,
            err_msg="scipy-NNLS projection must match sklearn's NMF.transform")
