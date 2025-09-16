#!/usr/bin/env python3
"""
Differential End Motif Analysis
===============================

This script performs differential analysis of fragment end motifs between two
groups of cfDNA samples. It identifies motifs that are significantly
over-represented in one group versus another using statistical testing.

The analysis uses Wilcoxon rank-sum tests to compare motif abundances between
groups, followed by FDR (False Discovery Rate) multiple testing correction
to control for the large number of motifs tested.

Input format:
    JSON configuration file with two fields:
    - "group_a": List of CSV file paths for group A samples
    - "group_b": List of CSV file paths for group B samples

    CSV files must have end motif format with columns:
    - motif_5p, count_5p, motif_3p, count_3p

Output:
    - Differential analysis results as CSV
    - Volcano plot showing effect sizes and significance
    - Summary statistics and significant motifs

Usage:
    python differential_end_motifs.py -c config.json -o output_dir/

Statistical Methods:
    - Wilcoxon rank-sum test for non-parametric comparison
    - Benjamini-Hochberg FDR correction for multiple testing
    - Effect size calculation using rank-biserial correlation

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
import json
import logging
import os
import sys
import pyfraglib

import numpy as np
import numpy.typing as npt
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from typing import Final, NoReturn
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection  # type: ignore

version_string: Final[str] = \
    "differential_end_motifs v{} (running on Python v{})" \
    .format(pyfraglib.__version__, sys.version.split(" ")[0])
LOGGER_NAME: Final[str] = "differential_end_motifs"


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


def load_config(
    config_path: str, logger: logging.Logger
) -> dict[str, list[str]]:
    """
    Load JSON configuration file with group definitions.

    Args:
        config_path: Path to JSON configuration file
        logger: Logger instance for output

    Returns:
        Dictionary with 'group_a' and 'group_b' file lists

    Raises:
        SystemExit: If config file is invalid or missing required fields
    """
    if not os.path.exists(config_path):
        fail("config file '{}' does not exist".format(config_path), logger)

    try:
        with open(config_path, "r") as f:
            config: dict[str, object] = json.load(f)
        required_fields: list[str] = ["group_a", "group_b"]
        for field in required_fields:
            if field not in config:
                fail("missing required field '{}' in config".format(field),
                     logger)

        group_a: list[str] = config["group_a"]  # type: ignore
        group_b: list[str] = config["group_b"]  # type: ignore
        if not isinstance(group_a, list) or not isinstance(group_b, list):
            fail("group_a and group_b must be lists of file paths", logger)

        if len(group_a) == 0 or len(group_b) == 0:
            fail("both groups must contain at least one file", logger)

        all_files: list[str] = group_a + group_b
        for file_path in all_files:
            if not os.path.exists(file_path):
                fail("file '{}' does not exist".format(file_path), logger)

        logger.info("found files: {} in group A, {} in group B".format(
            len(group_a), len(group_b)))
        return {"group_a": group_a, "group_b": group_b}

    except json.JSONDecodeError as e:
        fail("invalid JSON in config file: {}".format(e), logger)
    except Exception as e:
        fail("error loading config file: {}".format(e), logger)


def load_end_motif_data(
    csv_files: list[str], group_name: str, logger: logging.Logger
) -> pd.DataFrame:
    """
    Load end motif data from CSV files into a combined matrix.

    Args:
        csv_files: List of CSV file paths
        group_name: Name of the group (for logging)
        logger: Logger instance for output

    Returns:
        DataFrame with samples as rows and motifs as columns
        Contains both 5p and 3p motif counts

    Raises:
        SystemExit: If CSV files are malformed or incompatible
    """
    logger.info("loading {} files for group '{}'".format(
        len(csv_files), group_name))

    sample_data: dict[str, pd.DataFrame] = {}
    all_motifs_5p: set[str] = set()
    all_motifs_3p: set[str] = set()

    for csv_file in csv_files:
        sample_name: str = Path(csv_file).stem
        logger.debug("loading data from '{}'".format(csv_file))

        try:
            df: pd.DataFrame = pd.read_csv(csv_file)
            required_cols: list[str] = [
                "motif_5p", "count_5p", "motif_3p", "count_3p"
            ]
            if not all(col in df.columns for col in required_cols):
                fail("CSV file '{}' missing required columns: {}".format(
                    csv_file, required_cols), logger)

            df["count_5p"] = df["count_5p"].astype(int)  # type: ignore
            df["count_3p"] = df["count_3p"].astype(int)  # type: ignore
            df = df[
                (df["count_5p"] > 0) | (df["count_3p"] > 0)  # type: ignore
            ]

            if df.empty:
                fail("no valid data found in CSV file '{}'".format(
                    csv_file), logger)

            sample_data[sample_name] = df
            all_motifs_5p.update(df["motif_5p"].tolist())  # type: ignore
            all_motifs_3p.update(df["motif_3p"].tolist())  # type: ignore

        except Exception as e:
            fail("error reading CSV file '{}': {}".format(csv_file, e), logger)

    all_motifs_5p_sorted: list[str] = sorted(all_motifs_5p)
    all_motifs_3p_sorted: list[str] = sorted(all_motifs_3p)
    logger.info("found {} unique 5p motifs, {} unique 3p motifs".format(
        len(all_motifs_5p_sorted), len(all_motifs_3p_sorted)))

    sample_names: list[str] = list(sample_data.keys())
    n_samples: int = len(sample_names)
    n_motifs_total: int = len(all_motifs_5p_sorted) + len(all_motifs_3p_sorted)
    column_names: list[str] = []
    column_names.extend(
        ["motif_5p_{}".format(m) for m in all_motifs_5p_sorted]
    )
    column_names.extend(
        ["motif_3p_{}".format(m) for m in all_motifs_3p_sorted]
    )

    data_matrix: npt.NDArray[np.int64] = \
        np.zeros((n_samples, n_motifs_total), dtype=np.int64)

    for i, sample_name in enumerate(sample_names):
        df = sample_data[sample_name]
        for _, row in df.iterrows():  # type: ignore
            motif_5p: str = row["motif_5p"]  # type: ignore
            count_5p: int = row["count_5p"]  # type: ignore
            if motif_5p in all_motifs_5p_sorted and count_5p > 0:
                col_idx: int = all_motifs_5p_sorted.index(motif_5p)
                data_matrix[i, col_idx] = count_5p

            motif_3p: str = row["motif_3p"]  # type: ignore
            count_3p: int = row["count_3p"]  # type: ignore
            if motif_3p in all_motifs_3p_sorted and count_3p > 0:
                col_idx = len(all_motifs_5p_sorted) + \
                    all_motifs_3p_sorted.index(motif_3p)
                data_matrix[i, col_idx] = count_3p

    data_df: pd.DataFrame = pd.DataFrame(
        data_matrix, index=sample_names, columns=column_names
    )
    logger.debug("created data matrix: {} samples × {} motifs".format(
        data_df.shape[0], data_df.shape[1]))

    return data_df


def perform_differential_analysis(
    group_a_data: pd.DataFrame, group_b_data: pd.DataFrame,
    logger: logging.Logger
) -> pd.DataFrame:
    """
    Perform differential analysis between two groups using Wilcoxon tests.

    Args:
        group_a_data: DataFrame with group A samples and motif counts
        group_b_data: DataFrame with group B samples and motif counts
        logger: Logger instance for output

    Returns:
        DataFrame with differential analysis results for each motif

    Statistical Methods:
        - Wilcoxon rank-sum test (Mann-Whitney U) for non-parametric comparison
        - Benjamini-Hochberg FDR correction for multiple testing
        - Effect size using rank-biserial correlation
    """
    common_motifs: list[str] = \
        list(set(group_a_data.columns) & set(group_b_data.columns))
    common_motifs.sort()

    if len(common_motifs) == 0:
        fail("no common motifs found between groups", logger)
    logger.info("analyzing {} common motifs".format(len(common_motifs)))

    results: list[dict[str, object]] = []
    for motif in common_motifs:
        group_a_values: npt.NDArray[np.int64] = \
            group_a_data[motif].values  # type: ignore
        group_b_values: npt.NDArray[np.int64] = \
            group_b_data[motif].values  # type: ignore

        mean_a: float = float(np.mean(group_a_values))  # type: ignore
        mean_b: float = float(np.mean(group_b_values))  # type: ignore
        std_a: float = float(np.std(group_a_values))  # type: ignore
        std_b: float = float(np.std(group_b_values))  # type: ignore
        log_fold_change: float = np.log2((mean_b + 1) / (mean_a + 1))

        try:
            statistic, p_value = stats.mannwhitneyu(
                group_a_values, group_b_values, alternative="two-sided"
            )
            n_a: int = len(group_a_values)
            n_b: int = len(group_b_values)
            u_statistic: float = float(statistic)
            effect_size: float = 1 - (2 * u_statistic) / (n_a * n_b)
        except ValueError:
            p_value = 1.0
            effect_size = 0.0

        results.append({
            "motif": motif,
            "mean_group_a": mean_a,
            "mean_group_b": mean_b,
            "std_group_a": std_a,
            "std_group_b": std_b,
            "log_fold_change": log_fold_change,
            "p_value": p_value,
            "effect_size": effect_size,
        })

    results_df: pd.DataFrame = pd.DataFrame(results)
    rejected, p_adjusted = \
        fdrcorrection(results_df["p_value"].values, alpha=0.05)  # type: ignore
    results_df["p_adjusted"] = p_adjusted  # type: ignore
    results_df["significant"] = rejected  # type: ignore
    results_df = results_df.sort_values("p_value")
    n_significant: int = int(results_df["significant"].sum())  # type: ignore

    logger.info("found {} significant motifs after FDR correction".format(
        n_significant))
    return results_df


def create_differential_plots(
    results_df: pd.DataFrame, out_dir: str, logger: logging.Logger
) -> None:
    """
    Create visualization plots for differential analysis results.

    Args:
        results_df: DataFrame with differential analysis results
        out_dir: Output directory for plots
        logger: Logger instance for output
    """
    neg_log10_p_adj: npt.NDArray[np.float64] = \
        -np.log10(results_df["p_adjusted"].values + 1e-100)  # type: ignore
    colors: list[str] = [  # type: ignore
        "red" if sig else "gray"  # type: ignore
        for sig in results_df["significant"]  # type: ignore
    ]

    fig1, ax1 = plt.subplots(figsize=(10, 8))
    ax1.scatter(
        results_df["log_fold_change"], neg_log10_p_adj,  # type: ignore
        c=colors, alpha=0.7, s=30
    )

    ax1.axhline(
        y=-np.log10(0.05), color="black",  # type: ignore
        linestyle="--", alpha=0.5, label="FDR = 0.05"
    )
    ax1.axvline(x=0, color="black", linestyle="-", alpha=0.3)
    ax1.set_xlabel("Log2 Fold Change (Group B / Group A)")
    ax1.set_ylabel("-log10(FDR-adjusted p-value)")
    ax1.set_title("Volcano Plot: Differential End Motif Analysis")
    ax1.legend(["FDR = 0.05", "FDR < 0.05", "FDR ≥ 0.05"])  # type: ignore
    ax1.grid(True, alpha=0.3)

    top_significant: pd.DataFrame = \
        results_df[results_df["significant"]].head(15)  # type: ignore
    for _, row in top_significant.iterrows():  # type: ignore
        if row["p_adjusted"] < 0.01:  # type: ignore
            motif_name: str = (
                row["motif"].replace("motif_", "")  # type: ignore
                            .replace("_", " ")
            )
            ax1.annotate(
                motif_name,
                (row["log_fold_change"],  # type: ignore
                 -np.log10(row["p_adjusted"] + 1e-100)),  # type: ignore
                xytext=(5, 5), textcoords="offset points",
                fontsize=8, alpha=0.8
            )
    plt.tight_layout()

    volcano_path: str = os.path.join(out_dir, "differential_volcano_plot.png")
    fig1.savefig(volcano_path, dpi=300, bbox_inches="tight")
    logger.info("saved volcano plot to '{}'".format(volcano_path))
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    sig_effects: npt.NDArray[np.float64] = \
        results_df[
            results_df["significant"]  # type: ignore
        ]["effect_size"].values
    nonsig_effects: npt.NDArray[np.float64] = \
        results_df[
            ~results_df["significant"]  # type: ignore
        ]["effect_size"].values

    ax2.hist(nonsig_effects, bins=30, alpha=0.7,
             label="Non-significant", color="gray")
    ax2.hist(sig_effects, bins=30, alpha=0.7, label="Significant", color="red")

    ax2.set_xlabel("Effect Size (Rank-biserial correlation)")
    ax2.set_ylabel("Number of Motifs")
    ax2.set_title("Distribution of Effect Sizes")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    effect_path: str = os.path.join(out_dir, "differential_effect_sizes.png")
    fig2.savefig(effect_path, dpi=300, bbox_inches="tight")
    logger.info("saved effect size plot to '{}'".format(effect_path))
    plt.close(fig2)

    if len(top_significant) > 0:
        fig3, ax3 = plt.subplots(figsize=(12, 8))
        top_20: pd.DataFrame = \
            results_df[results_df["significant"]].head(20)  # type: ignore
        motif_labels: list[str] = [  # type: ignore
            m.replace("motif_", "").replace("_", " ")  # type: ignore
            for m in top_20["motif"].tolist()  # type: ignore
        ]
        bars = ax3.barh(
            range(len(top_20)),
            top_20["log_fold_change"].values  # type: ignore
        )

        for i, bar in enumerate(bars):  # type: ignore
            if top_20.iloc[i]["log_fold_change"] > 0:  # type: ignore
                bar.set_color("red")  # type: ignore
            else:
                bar.set_color("blue")  # type: ignore

        ax3.set_yticks(range(len(top_20)))  # type: ignore
        ax3.set_yticklabels(motif_labels)  # type: ignore
        ax3.set_xlabel("Log2 Fold Change")
        ax3.set_title("Top 20 Significant Differential Motifs")
        ax3.axvline(x=0, color="black", linestyle="-", alpha=0.5)
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        top_path: str = os.path.join(out_dir, "differential_top_motifs.png")
        fig3.savefig(top_path, dpi=300, bbox_inches="tight")
        logger.info("saved top motifs plot to '{}'".format(top_path))
        plt.close(fig3)


def create_argparser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser.

    Returns:
        Configured ArgumentParser instance
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="differential_end_motifs",
        description="Perform differential analysis of end motif abundances.",
        epilog="{}. Licensed under GPLv3.".format(version_string))

    parser.add_argument(
        "-c", "--config", type=str, dest="config_file", required=True,
        help="JSON configuration file with group_a and group_b file lists")

    parser.add_argument(
        "-o", "--out-dir", type=str, dest="out_dir",
        default="differential_output",
        help="Output directory for results and plots (default: "
        "differential_output)")

    parser.add_argument(
        "-v", "--verbose", action="store_true", dest="is_verbose",
        default=False, help="Enable verbose logging")

    return parser


def main() -> None:
    """
    Main function for the differential end motif analysis script.
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
    logger.info("starting differential end motif analysis")
    logger.info(version_string)

    out_dir: str = args.out_dir
    if not os.path.exists(out_dir):
        logger.info("creating output directory '{}'".format(out_dir))
        os.makedirs(out_dir, exist_ok=True)

    config_file: str = args.config_file
    config: dict[str, list[str]] = load_config(config_file, logger)
    group_a_data: pd.DataFrame = load_end_motif_data(
        config["group_a"], "A", logger
    )
    group_b_data: pd.DataFrame = load_end_motif_data(
        config["group_b"], "B", logger
    )

    results_df: pd.DataFrame = perform_differential_analysis(
        group_a_data, group_b_data, logger
    )
    results_path: str = os.path.join(out_dir, "differential_results.csv")
    results_df.to_csv(results_path, index=False)
    create_differential_plots(results_df, out_dir, logger)


if __name__ == "__main__":
    main()
