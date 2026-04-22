#!/usr/bin/env python3
"""
Non-negative Matrix Factorization of Fragment Length Distributions
==================================================================

This script performs non-negative matrix factorization (NMF) on fragment length
distributions from cfDNA samples. NMF can identify underlying patterns or
signatures in those length distributions.

NMF decomposes the fragment length matrix into basis components (signatures)
and mixing coefficients (weights), which can reveal tissue-specific or
pathological patterns in cfDNA fragmentation.

Input format:
- Text file with one fragment length CSV file path per line
- CSV files must have columns: fragment_length, count

Output:
- NMF basis matrix (signatures) as CSV
- NMF coefficients matrix (weights) as CSV
- Visualization plots showing signatures and sample compositions
- Reconstruction error analysis

Usage:
    # Fit mode:
    python nmf_fragment_lengths.py -f file_list.txt -n 3 -o output_dir/

    # Projection mode (onto an existing nmf_signatures.csv):
    python nmf_fragment_lengths.py -f new_samples.txt \
        -s output_dir/nmf_signatures.csv -o projected_dir/

Parameters:
    -n, --n-components: Number of NMF components (signatures) to extract
    -o, --out-dir: Output directory for results and plots
    -f, --file-list: Text file containing CSV file paths
    -s, --signatures: Optional path to an existing nmf_signatures.csv. When
        set, the script skips refitting and projects the input samples onto
        the loaded signatures basis (equivalent to a fixed-H NMF transform).

License
-------
Copyright (C) 2026 Daniel Schütte, daniel.schuette@iccb-cologne.org

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details. You should have received a copy of the GNU General Public
License along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
import argparse
import logging
import os
import sys
import pyfraglib

import numpy as np
import numpy.typing as npt
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from typing import Final, NoReturn
from scipy.optimize import nnls
from sklearn.decomposition import NMF  # type: ignore
from sklearn.metrics import mean_squared_error  # type: ignore

version_string: Final[str] = \
    "nmf_fragment_lengths v{} (running on Python v{})" \
    .format(pyfraglib.__version__, sys.version.split(" ")[0])
LOGGER_NAME: Final[str] = "nmf_fragment_lengths"


def fail(msg: str, logger: logging.Logger) -> NoReturn:
    """
    Log a fatal error message and exit the program.

    Args:
        msg: Error message to log
        logger: Logger instance to use for logging

    Raises:
        SystemExit: Always exits with code 1
    """
    logger.fatal(msg)
    sys.exit(1)


def read_file_list(file_list_path: str, logger: logging.Logger) -> list[str]:
    """
    Read list of CSV file paths from text file.

    Args:
        file_list_path: Path to text file containing CSV file paths
        logger: Logger instance for output

    Returns:
        List of CSV file paths

    Raises:
        SystemExit: If file doesn't exist or is empty
    """
    if not os.path.exists(file_list_path):
        fail("file list '{}' does not exist".format(file_list_path), logger)

    csv_files: list[str] = []
    with open(file_list_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if not os.path.exists(line):
                fail("CSV file '{}' (line {}) does not exist".format(
                    line, line_num), logger)
            csv_files.append(line)

    if not csv_files:
        fail("no valid CSV files found in file list", logger)
    logger.info("found {} CSV files to process".format(len(csv_files)))

    return csv_files


def load_fragment_length_data(
    csv_files: list[str], logger: logging.Logger
) -> tuple[pd.DataFrame, list[str]]:
    """
    Load fragment length data from CSV files into a matrix.

    Args:
        csv_files: List of CSV file paths
        logger: Logger instance for output

    Returns:
        Tuple of (data matrix, sample names)
        Data matrix has samples as rows and fragment lengths as columns

    Raises:
        SystemExit: If CSV files are malformed or incompatible
    """
    logger.info("loading fragment length data from {} files".format(
        len(csv_files)))

    sample_data: dict[str, pd.DataFrame] = {}
    all_lengths: set[int] = set()

    for csv_file in csv_files:
        sample_name: str = Path(csv_file).stem
        logger.debug("loading data from '{}'".format(csv_file))

        try:
            df: pd.DataFrame = pd.read_csv(csv_file)

            required_cols: list[str] = ["fragment_length", "count"]
            if not all(col in df.columns for col in required_cols):
                fail("CSV file '{}' missing required columns: {}".format(
                    csv_file, required_cols), logger)

            df["fragment_length"] = df["fragment_length"].astype(int)
            df["count"] = df["count"].astype(int)
            df = df[(df["count"] > 0) & (df["fragment_length"] > 0)]

            if df.empty:
                fail("no valid data found in CSV file '{}'".format(
                    csv_file), logger)

            sample_data[sample_name] = df
            all_lengths.update(df["fragment_length"].tolist())

        except Exception as e:
            fail("error reading CSV file '{}': {}".format(csv_file, e), logger)

    min_length: int = min(all_lengths)
    max_length: int = max(all_lengths)
    length_range: range = range(min_length, max_length + 1)
    logger.info("fragment length range: {}-{} bp ({} total lengths)".format(
        min_length, max_length, len(length_range)))

    sample_names: list[str] = list(sample_data.keys())
    data_matrix: npt.NDArray[np.float64] = \
        np.zeros((len(sample_names), len(length_range)))

    for i, sample_name in enumerate(sample_names):
        df = sample_data[sample_name]
        row: pd.Series[int]
        for _, row in df.iterrows():
            length_idx: int = row["fragment_length"] - min_length
            data_matrix[i, length_idx] = row["count"]

    column_names: list[str] = [str(length) for length in length_range]
    data_df: pd.DataFrame = pd.DataFrame(
        data_matrix, index=sample_names, columns=column_names
    )
    logger.info("created data matrix: {} samples × {} fragment lengths".format(
        data_df.shape[0], data_df.shape[1]))

    return data_df, sample_names


def load_signatures(
    signatures_path: str, logger: logging.Logger,
) -> tuple[npt.NDArray[np.float64], list[str]]:
    """
    Load a previously-fit NMF signatures matrix (``H``).

    Args:
        signatures_path: Path to an ``nmf_signatures.csv`` produced by a
            previous run. Rows index components, columns index fragment
            lengths (integer bp labels).
        logger: Logger instance for output.

    Returns:
        Tuple of the ``H`` matrix (components × lengths) and the list of
        fragment-length column labels.

    Raises:
        SystemExit: If the file does not exist or is malformed.
    """
    if not os.path.exists(signatures_path):
        fail("signatures file '{}' does not exist".format(
            signatures_path), logger)

    try:
        df: pd.DataFrame = pd.read_csv(signatures_path, index_col=0)
    except Exception as e:
        fail("error reading signatures file '{}': {}".format(
            signatures_path, e), logger)

    lengths: list[str] = [str(c) for c in df.columns]
    H: npt.NDArray[np.float64] = df.to_numpy(dtype=np.float64)
    if H.ndim != 2 or H.shape[0] < 1 or H.shape[1] < 1:
        fail("signatures file '{}' has unexpected shape {}".format(
            signatures_path, H.shape), logger)

    logger.info(
        "loaded signatures: {} components × {} fragment lengths".format(
            H.shape[0], H.shape[1]
        )
    )
    return H, lengths


def align_to_signatures(
    data_matrix: pd.DataFrame, signature_lengths: list[str],
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Align a sample x fragment-length matrix to the column layout of an
    existing signatures matrix, filling missing lengths with zero and
    dropping any extra lengths.

    Args:
        data_matrix: Sample × fragment-length matrix from
            :func:`load_fragment_length_data`.
        signature_lengths: Fragment-length column labels from the signatures
            file.
        logger: Logger instance for output.

    Returns:
        Aligned matrix with the same columns (and order) as the signatures
        file.
    """
    aligned: pd.DataFrame = data_matrix.reindex(
        columns=signature_lengths, fill_value=0.0
    )
    n_missing: int = int((aligned.sum(axis=0) == 0).sum())
    n_dropped: int = int(len(data_matrix.columns) - len(signature_lengths))
    if n_missing:
        logger.warning(
            "{} fragment-length columns from the signatures basis have no "
            "counts in the projection samples (filled with 0)".format(
                n_missing))
    if n_dropped > 0:
        logger.warning(
            "{} fragment-length columns present in the projection samples "
            "are absent from the signatures basis and were dropped".format(
                n_dropped))
    return aligned


def project_onto_signatures(
    data_matrix: pd.DataFrame, H: npt.NDArray[np.float64],
    logger: logging.Logger,
) -> tuple[npt.NDArray[np.float64], float]:
    """
    Project row-normalised samples onto a fixed signatures matrix.

    Solves, for each sample :math:`x_i`,

    .. math::

        w_i = \\arg\\min_{w \\ge 0} \\| x_i - w\\, H \\|_2^2

    using :func:`scipy.optimize.nnls`. This is the non-negative least
    squares dual of NMF with :math:`H` held fixed, equivalent to what
    :meth:`sklearn.decomposition.NMF.transform` does internally when the
    coefficient matrix is frozen, and avoids the need for private sklearn
    APIs.

    Args:
        data_matrix: Sample × fragment-length matrix aligned to the
            signatures column layout.
        H: Fixed signatures matrix (components × lengths).
        logger: Logger instance for output.

    Returns:
        Tuple of (``W``, reconstruction error). ``W`` has shape
        ``(n_samples, n_components)``.
    """
    raw: npt.NDArray[np.float64] = data_matrix.values.astype(np.float64)
    row_sums: npt.NDArray[np.float64] = raw.sum(axis=1, keepdims=True)
    X_scaled: npt.NDArray[np.float64] = raw / np.where(
        row_sums > 0, row_sums, 1.0
    )

    n_samples: int = X_scaled.shape[0]
    n_components: int = H.shape[0]
    W: npt.NDArray[np.float64] = np.zeros((n_samples, n_components))
    H_T: npt.NDArray[np.float64] = H.T
    for i in range(n_samples):
        w, _ = nnls(H_T, X_scaled[i])
        W[i] = w

    X_reconstructed: npt.NDArray[np.float64] = np.dot(W, H)
    error: float = mean_squared_error(X_scaled, X_reconstructed)

    logger.info("projection reconstruction error: {:.6f}".format(error))
    return W, error


def perform_nmf_analysis(
    data_matrix: pd.DataFrame, n_components: int, logger: logging.Logger
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], float]:
    """
    Perform non-negative matrix factorization on fragment length data.

    Args:
        data_matrix: Sample×length data matrix
        n_components: Number of NMF components to extract
        logger: Logger instance for output

    Returns:
        Tuple of (W, H, reconstruction error). W is the raw (unnormalized)
        sample loading matrix from NMF. To obtain true mixture fractions,
        weight each column of W by the corresponding H row sum before
        row-normalizing (H-weighted normalization), rather than plain
        row-normalization of W, which ignores differences in signature
        magnitudes.
    """
    logger.info("performing NMF with {} components".format(n_components))

    raw: npt.NDArray[np.float64] = data_matrix.values.astype(np.float64)
    row_sums: npt.NDArray[np.float64] = raw.sum(axis=1, keepdims=True)
    X_scaled: npt.NDArray[np.float64] = raw / np.where(
        row_sums > 0, row_sums, 1.0
    )

    # @NOTE(ds): We perform NMF with multiple random initializations for
    # better stability.
    best_error: float = float("inf")
    best_W: npt.NDArray[np.float64] | None = None
    best_H: npt.NDArray[np.float64] | None = None

    n_runs: int = 10
    logger.info("running NMF {} times with different initializations".format(
        n_runs))

    for run in range(n_runs):
        nmf_model: NMF = NMF(  # type: ignore
            n_components=n_components,
            init="random",
            random_state=run,
            max_iter=1000,
        )

        W: npt.NDArray[np.float64] = nmf_model.fit_transform(X_scaled)
        H: npt.NDArray[np.float64] = nmf_model.components_
        X_reconstructed: npt.NDArray[np.float64] = np.dot(W, H)
        error: float = mean_squared_error(X_scaled, X_reconstructed)

        if error < best_error:
            best_error = error
            best_W = W
            best_H = H

        logger.debug("NMF run {}: reconstruction error = {:.6f}".format(
            run + 1, error))

    logger.info("best reconstruction error: {:.6f}".format(best_error))
    if best_W is None or best_H is None:
        fail("NMF failed to converge", logger)

    return best_W, best_H, best_error


def save_nmf_results(
    W: npt.NDArray[np.float64], H: npt.NDArray[np.float64],
    sample_names: list[str], fragment_lengths: list[str], out_dir: str,
    logger: logging.Logger
) -> None:
    """
    Save NMF results to a CSV files.

    Args:
        W: Basis matrix (samples × components)
        H: Coefficient matrix (components × fragment lengths)
        sample_names: List of sample names
        fragment_lengths: List of fragment length strings
        out_dir: Output directory
        logger: Logger instance for output
    """
    component_names: list[str] = [
        "Component_{}".format(i+1) for i in range(W.shape[1])
    ]
    W_df: pd.DataFrame = pd.DataFrame(
        W, index=sample_names, columns=component_names
    )
    w_path: str = os.path.join(out_dir, "nmf_mixing_coefficients.csv")
    W_df.to_csv(w_path)
    logger.info("saved mixing coefficients to '{}'".format(w_path))

    H_df: pd.DataFrame = pd.DataFrame(
        H, index=component_names, columns=fragment_lengths
    )
    h_path: str = os.path.join(out_dir, "nmf_signatures.csv")
    H_df.to_csv(h_path)
    logger.info("saved signatures to '{}'".format(h_path))


def create_nmf_visualizations(
    W: npt.NDArray[np.float64], H: npt.NDArray[np.float64],
    sample_names: list[str], fragment_lengths: list[str],
    reconstruction_error: float, out_dir: str, logger: logging.Logger
) -> None:
    """
    Create comprehensive visualizations of NMF results.

    Args:
        W: Basis matrix (samples × components)
        H: Coefficient matrix (components × fragment lengths)
        sample_names: List of sample names
        fragment_lengths: List of fragment length strings
        reconstruction_error: NMF reconstruction error
        out_dir: Output directory
        logger: Logger instance for output
    """
    n_components: int = W.shape[1]
    component_names: list[str] = [
        "Component {}".format(i+1) for i in range(n_components)
    ]
    fig1, axes1 = plt.subplots(
        n_components, 1, figsize=(12, 3*n_components), sharex=True
    )
    if n_components == 1:
        axes1 = [axes1]

    lengths_numeric: list[int] = [int(x) for x in fragment_lengths]
    ax: plt.Axes
    i: int
    for i, ax in enumerate(axes1):
        ax.plot(
            lengths_numeric, H[i, :], linewidth=2, label=component_names[i]
        )
        ax.set_ylabel("Signature Weight")
        ax.set_title("NMF Component {} Signature".format(i+1))
        ax.grid(True, alpha=0.3)
        ax.legend()

    axes1[-1].set_xlabel("Fragment Length (bp)")
    plt.suptitle("NMF Fragment Length Signatures\n"
                 "Reconstruction Error: {:.6f}".format(reconstruction_error))
    plt.tight_layout()

    sig_path: str = os.path.join(out_dir, "nmf_signatures.png")
    fig1.savefig(sig_path, dpi=300, bbox_inches="tight")
    logger.info("saved signatures plot to '{}'".format(sig_path))
    plt.close(fig1)

    ax2: plt.Axes
    fig2, ax2 = plt.subplots(figsize=(max(8, len(sample_names)*0.5), 6))
    W_df: pd.DataFrame = pd.DataFrame(
        W, index=sample_names, columns=component_names
    )
    # H-weighted normalisation: scale each W column by the corresponding
    # H row sum before row-normalising, so that signature magnitude differences
    # are accounted for and the result represents true mixture fractions.
    H_sums: npt.NDArray[np.float64] = H.sum(axis=1)
    W_weighted: pd.DataFrame = W_df.multiply(H_sums, axis="columns")
    W_norm: pd.DataFrame = W_weighted.div(W_weighted.sum(axis=1), axis=0)
    sns.heatmap(
        W_norm.T, annot=True, fmt=".3f", cmap="viridis",
        ax=ax2, cbar_kws={"label": "Mixture fraction"}
    )
    ax2.set_title("Sample Composition by NMF Components")
    ax2.set_xlabel("Samples")
    ax2.set_ylabel("NMF Components")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    comp_path: str = os.path.join(out_dir, "nmf_sample_composition.png")
    fig2.savefig(comp_path, dpi=300, bbox_inches="tight")
    logger.info("saved composition plot to '{}'".format(comp_path))
    plt.close(fig2)

    ax3: plt.Axes
    fig3, ax3 = plt.subplots(figsize=(max(8, len(sample_names)*0.5), 6))
    x_pos: npt.NDArray[np.int64] = np.arange(len(sample_names))
    width: float = 0.8 / n_components

    for i in range(n_components):
        ax3.bar(
            x_pos + i*width, W[:, i], width,
            label=component_names[i], alpha=0.8
        )

    ax3.set_xlabel("Samples")
    ax3.set_ylabel("Component Weight")
    ax3.set_title("NMF Component Weights by Sample")
    ax3.set_xticks(x_pos + width*(n_components-1)/2)
    ax3.set_xticklabels(sample_names, rotation=45, ha="right")  # type: ignore
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    plt.tight_layout()

    weights_path: str = os.path.join(out_dir, "nmf_component_weights.png")
    fig3.savefig(weights_path, dpi=300, bbox_inches="tight")
    logger.info("saved weights plot to '{}'".format(weights_path))
    plt.close(fig3)


def create_argparser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser.

    Returns:
        Configured ArgumentParser instance
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="nmf_fragment_lengths",
        description="Perform NMF analysis on fragment length distributions.",
        epilog="{}. Licensed under GPLv3.".format(version_string))

    parser.add_argument(
        "-f", "--file-list", type=str, dest="file_list", required=True,
        help="Text file containing fragment length CSV file paths "
        "(one per line)")

    parser.add_argument(
        "-n", "--n-components", type=int, dest="n_components", default=3,
        help="Number of NMF components to extract (default: 3). Ignored "
        "when --signatures is given; the component count is then taken "
        "from the loaded signatures file.")

    parser.add_argument(
        "-o", "--out-dir", type=str, dest="out_dir", default="nmf_output",
        help="Output directory for results and plots (default: nmf_output)")

    parser.add_argument(
        "-s", "--signatures", type=str, dest="signatures_file",
        default=None, help="Optional path to an existing nmf_signatures.csv "
        "from a previous run. When set, the script skips refitting and "
        "instead projects the input samples onto the loaded signatures "
        "basis (equivalent to sklearn's NMF.transform with a fixed H). "
        "Fragment-length columns in the projection samples are aligned to "
        "the basis layout before projection.")

    parser.add_argument(
        "-v", "--verbose", action="store_true", dest="is_verbose",
        default=False, help="Enable verbose logging")

    return parser


def main() -> None:
    """
    Main function for the NMF fragment length analysis script.
    """
    logger: logging.Logger = logging.getLogger(LOGGER_NAME)
    parser: argparse.ArgumentParser = create_argparser()
    args: argparse.Namespace = parser.parse_args()

    level: int
    is_verbose: Final[bool] = args.is_verbose
    if is_verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logger.setLevel(level)

    n_components: int = args.n_components
    if n_components < 1:
        fail("number of components must be >= 1", logger)

    out_dir: str = args.out_dir
    if not os.path.exists(out_dir):
        logger.info("creating output directory '{}'".format(out_dir))
        os.makedirs(out_dir, exist_ok=True)

    file_list: str = args.file_list
    csv_files: list[str] = read_file_list(file_list, logger)
    data_matrix, sample_names = load_fragment_length_data(csv_files, logger)

    signatures_file: str | None = args.signatures_file
    if signatures_file is not None:
        H, signature_lengths = load_signatures(signatures_file, logger)
        data_matrix = align_to_signatures(
            data_matrix, signature_lengths, logger
        )
        W, reconstruction_error = project_onto_signatures(
            data_matrix, H, logger
        )
        fragment_lengths: list[str] = signature_lengths
    else:
        n_components: int = int(args.n_components)  # type: ignore
        W, H, reconstruction_error = perform_nmf_analysis(
            data_matrix, n_components, logger
        )
        fragment_lengths = data_matrix.columns.tolist()

    save_nmf_results(W, H, sample_names, fragment_lengths, out_dir, logger)
    create_nmf_visualizations(
        W, H, sample_names, fragment_lengths, reconstruction_error,
        out_dir, logger)


if __name__ == "__main__":
    main()
