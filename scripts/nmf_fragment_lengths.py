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
    python nmf_fragment_lengths.py -f file_list.txt -n 3 -o output_dir/

Parameters:
    -n, --n-components: Number of NMF components (signatures) to extract
    -o, --out-dir: Output directory for results and plots
    -f, --file-list: Text file containing CSV file paths

License
-------
Copyright (C) 2025 Daniel Schütte, daniel.schuette@iccb-cologne.org

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
        Tuple of (basis matrix W, coefficient matrix H, reconstruction error)
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

    # @NOTE(ds): We Aormalize W matrix so each row (sample) sums to 1.
    # This makes coefficients represent true proportions of each component.
    W_row_sums: npt.NDArray[np.float64] = best_W.sum(axis=1, keepdims=True)
    W_normalized: npt.NDArray[np.float64] = best_W / W_row_sums
    return W_normalized, best_H, best_error


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
    W_norm: pd.DataFrame = W_df.div(W_df.sum(axis=1), axis=0)
    sns.heatmap(
        W_norm.T, annot=True, fmt=".3f", cmap="viridis",
        ax=ax2, cbar_kws={"label": "Normalized Weight"}
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
        help="Number of NMF components to extract (default: 3)")

    parser.add_argument(
        "-o", "--out-dir", type=str, dest="out_dir", default="nmf_output",
        help="Output directory for results and plots (default: nmf_output)")

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

    n_components: int = int(args.n_components)  # type: ignore
    W, H, reconstruction_error = perform_nmf_analysis(
        data_matrix, n_components, logger
    )

    fragment_lengths: list[str] = data_matrix.columns.tolist()
    save_nmf_results(W, H, sample_names, fragment_lengths, out_dir, logger)
    create_nmf_visualizations(
        W, H, sample_names, fragment_lengths, reconstruction_error,
        out_dir, logger)


if __name__ == "__main__":
    main()
