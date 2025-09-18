#!/usr/bin/env python3
"""
Nuclease Parameter Estimation Using ABC-SMC
===========================================

This script uses approximate Bayesian computation with sequential Monte-Carlo
to learn the nuclease parameters used by ``pyfraglib``s simulator from data.
It is mostly a proof of concept.

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
import pyabc  # type: ignore
import pyfraglib

import numpy as np
import pandas as pd

from typing import Final, NoReturn
from functools import partial
from pyfraglib import FragmentSimulator, NucleaseProfile, FragmentList
from pyfraglib.simulator.fragment_simulator import SequenceContextGenerator
from pyfraglib.core import detect_cpus

version_string: Final[str] = \
    "learn_nuclease_params v{} (running on Python v{})".format(
        pyfraglib.__version__, sys.version.split(" ")[0]
    )
LOGGER_NAME: Final[str] = "learn_nuclease_params"


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


def load_observed_motifs(
    csv_path: str, logger: logging.Logger
) -> dict[str, float]:
    """
    Load and normalize observed 5' end motif frequencies from stats CSV output.

    Args:
        csv_path: Path to CSV file from pyfrag.py stats command

    Returns:
        Dictionary mapping motif sequences to normalized frequencies
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Motif data file not found: {csv_path}")

    df = pd.read_csv(csv_path)
    required_cols = ["motif_5p", "count_5p"]
    for col in required_cols:
        if col not in df.columns:
            fail(f"Required column '{col}' not found in {csv_path}", logger)

    motif_counts: dict[str, int] = dict(
        zip(df["motif_5p"], df["count_5p"])  # type: ignore
    )
    total_count = sum(motif_counts.values())
    if total_count == 0:
        fail("No motif counts found in data", logger)

    motif_freqs: dict[str, float] = {
        motif: count / total_count
        for motif, count in motif_counts.items()
    }
    return motif_freqs


def simulate_motifs(
    params: dict[str, object], fasta_path: str
) -> dict[str, float]:
    """
    Simulate fragment 5' end motifs given nuclease parameters.

    Args:
        params: Dict containing nuclease activity and preference parameters

    Returns:
        Dict mapping motif sequences to simulated frequencies
    """
    dnase1_motif_preference: dict[str, float] = {
        "CG": max(params["dnase1_cg_pref"], 0.0),  # type: ignore
        "GC": max(params["dnase1_gc_pref"], 0.0),  # type: ignore
        "AT": max(params["dnase1_at_pref"], 0.0),  # type: ignore
        "TA": max(params["dnase1_ta_pref"], 0.0),  # type: ignore
        "AA": max(params["dnase1_aa_pref"], 0.0),  # type: ignore
        "TT": max(params["dnase1_tt_pref"], 0.0)   # type: ignore
    }

    dnase1l3_motif_preference: dict[str, float] = {
        "CC": max(params["dnase1l3_cc_pref"], 0.0),  # type: ignore
        "C": max(params["dnase1l3_c_pref"], 0.0),    # type: ignore
        "CG": max(params["dnase1l3_cg_pref"], 0.0),  # type: ignore
        "GC": max(params["dnase1l3_gc_pref"], 0.0),  # type: ignore
        "CT": max(params["dnase1l3_ct_pref"], 0.0),  # type: ignore
        "TC": max(params["dnase1l3_tc_pref"], 0.0),  # type: ignore
        "GG": max(params["dnase1l3_gg_pref"], 0.0),  # type: ignore
        "AT": max(params["dnase1l3_at_pref"], 0.0),  # type: ignore
        "TA": max(params["dnase1l3_ta_pref"], 0.0),  # type: ignore
        "A": max(params["dnase1l3_a_pref"], 0.0),    # type: ignore
        "T": max(params["dnase1l3_t_pref"], 0.0)     # type: ignore
    }

    dffb_motif_preference: dict[str, float] = {
        "A": max(params["dffb_a_pref"], 0.0),  # type: ignore
        "T": max(params["dffb_t_pref"], 0.0),  # type: ignore
        "C": max(params["dffb_c_pref"], 0.0),  # type: ignore
        "G": max(params["dffb_g_pref"], 0.0)   # type: ignore
    }

    dnase1_activity: float = \
        max(params["dnase1_activity"], 0.0)  # type: ignore
    dnase1l3_activity: float = \
        max(params["dnase1l3_activity"], 0.0)  # type: ignore
    dffb_activity: float = max(params["dffb_activity"], 0.0)  # type: ignore
    nuclease_profile = NucleaseProfile(
        dnase1_activity=dnase1_activity,
        dnase1l3_activity=dnase1l3_activity,
        dffb_activity=dffb_activity,
        dnase1_motif_preference=dnase1_motif_preference,
        dnase1l3_motif_preference=dnase1l3_motif_preference,
        dffb_motif_preference=dffb_motif_preference
    )

    seq_gen: SequenceContextGenerator = SequenceContextGenerator(fasta_path)
    try:
        simulator: FragmentSimulator = FragmentSimulator(seq_gen)
        fragments: FragmentList = simulator.simulate_fragments(
            chrom="1", start=100_000, end=200_000_000,
            num_fragments=50000, tissue_type="healthy",
            nuclease_profile=nuclease_profile
        )

        motifs_5p, _, num_frags = fragments.count_endmotifs(kmer_len=4)
        if num_frags == 0:
            raise ValueError("No 5' end motifs simulated")

        total_count: int = sum(motifs_5p.values())
        motif_freqs: dict[str, float] = {
            motif: count / total_count
            for motif, count in motifs_5p.items()
        }
        return motif_freqs
    finally:
        seq_gen.close()


def calculate_distance(
    simulated: dict[str, float], observed: dict[str, float]
) -> float:
    """
    Calculate L2 distance between simulated and observed motif frequencies.

    Args:
        simulated: Simulated motif frequencies
        observed: Observed motif frequencies from data

    Returns:
        Euclidean distance
    """
    all_motifs = set(simulated.keys()) | set(observed.keys())
    distance = 0.0
    for motif in all_motifs:
        sim_freq = simulated.get(motif, 0.0)
        obs_freq = observed.get(motif, 0.0)
        distance += (sim_freq - obs_freq) ** 2
    return float(np.sqrt(distance))  # type: ignore


def abc_model(
    params: dict[str, object], fasta_path: str
) -> dict[str, float]:
    """
    ABC model function that simulates motifs given parameters.
    """
    return simulate_motifs(params, fasta_path)


def abc_distance(
    simulated: dict[str, float], observed: dict[str, float]
) -> float:
    """
    ABC distance function that compares simulated to observed data.
    """
    return calculate_distance(simulated, observed)


def optimize(
    motif_data_path: str, fasta_path: str, output_dir: str,
    pop_size: int, logger: logging.Logger
) -> None:
    """
    Run ABC-SMC parameter estimation.

    Args:
        motif_data_path: Path to CSV file with observed motif data
        pop_size: Population size per iteration
        fasta_path: Path to reference FASTA file
        output_dir: Directory for output files
        logger: Logger instance
    """
    observed_motifs = load_observed_motifs(motif_data_path, logger)
    logger.info(f"Loaded {len(observed_motifs)} motifs from observed data")

    prior = pyabc.Distribution(  # type: ignore
        dnase1_activity=(
            pyabc.RV("lognorm", s=0.5, scale=1.0)  # type: ignore
        ),
        dnase1l3_activity=(
            pyabc.RV("lognorm", s=0.5, scale=1.2)  # type: ignore
        ),
        dffb_activity=(
            pyabc.RV("lognorm", s=0.8, scale=0.3)  # type: ignore
        ),

        dnase1_cg_pref=(
            pyabc.RV("lognorm", s=0.3, scale=0.8)  # type: ignore
        ),
        dnase1_gc_pref=(
            pyabc.RV("lognorm", s=0.3, scale=0.8)  # type: ignore
        ),
        dnase1_at_pref=(
            pyabc.RV("lognorm", s=0.25, scale=1.2)  # type: ignore
        ),
        dnase1_ta_pref=(
            pyabc.RV("lognorm", s=0.25, scale=1.2)  # type: ignore
        ),
        dnase1_aa_pref=(
            pyabc.RV("lognorm", s=0.2, scale=1.1)  # type: ignore
        ),
        dnase1_tt_pref=(
            pyabc.RV("lognorm", s=0.2, scale=1.1)  # type: ignore
        ),

        dnase1l3_cc_pref=(
            pyabc.RV("lognorm", s=0.4, scale=8.0)  # type: ignore
        ),
        dnase1l3_c_pref=(
            pyabc.RV("lognorm", s=0.35, scale=4.5)  # type: ignore
        ),
        dnase1l3_cg_pref=(
            pyabc.RV("lognorm", s=0.3, scale=3.5)  # type: ignore
        ),
        dnase1l3_gc_pref=(
            pyabc.RV("lognorm", s=0.3, scale=3.5)  # type: ignore
        ),
        dnase1l3_ct_pref=(
            pyabc.RV("lognorm", s=0.3, scale=2.8)  # type: ignore
        ),
        dnase1l3_tc_pref=(
            pyabc.RV("lognorm", s=0.3, scale=2.8)  # type: ignore
        ),
        dnase1l3_gg_pref=(
            pyabc.RV("lognorm", s=0.2, scale=1.0)  # type: ignore
        ),
        dnase1l3_at_pref=(
            pyabc.RV("lognorm", s=0.25, scale=0.6)  # type: ignore
        ),
        dnase1l3_ta_pref=(
            pyabc.RV("lognorm", s=0.25, scale=0.6)  # type: ignore
        ),
        dnase1l3_a_pref=(
            pyabc.RV("lognorm", s=0.2, scale=0.7)  # type: ignore
        ),
        dnase1l3_t_pref=(
            pyabc.RV("lognorm", s=0.2, scale=0.7)  # type: ignore
        ),

        dffb_a_pref=(
            pyabc.RV("lognorm", s=0.15, scale=1.1)  # type: ignore
        ),
        dffb_t_pref=(
            pyabc.RV("lognorm", s=0.15, scale=1.05)  # type: ignore
        ),
        dffb_c_pref=(
            pyabc.RV("lognorm", s=0.1, scale=0.95)  # type: ignore
        ),
        dffb_g_pref=(
            pyabc.RV("lognorm", s=0.1, scale=0.95)  # type: ignore
        )
    )

    n_cores = detect_cpus()
    logger.info(f"Using {n_cores} cores for ABC-SMC sampling")
    sampler = pyabc.sampler.MulticoreEvalParallelSampler(  # type: ignore
        n_procs=n_cores
    )

    # transition = pyabc.LocalTransition(  # type: ignore
    #     k_fraction=0.3, scaling=0.2
    # )
    transition = pyabc.MultivariateNormalTransition(  # type: ignore
        scaling=0.2
    )
    abc = pyabc.ABCSMC(  # type: ignore
        models=partial(abc_model, fasta_path=fasta_path),
        parameter_priors=prior,  # type: ignore
        distance_function=abc_distance,
        population_size=pop_size,
        transitions=transition,  # type: ignore
        sampler=sampler  # type: ignore
    )

    db_path = os.path.join(output_dir, "abc_results.db")
    logger.info(f"Initializing ABC database: {db_path}")
    history: pyabc.History = abc.new(  # type: ignore
        f"sqlite:///{db_path}", observed_motifs
    )

    logger.info("Running ABC-SMC inference...")
    history = abc.run(  # type: ignore
        minimum_epsilon=0.01,
        max_nr_populations=100,
        min_acceptance_rate=0.01
    )
    logger.info(f"Run results saved to ABC database: {db_path}")

    logger.info(f"Saving posterior estimates and plots to {output_dir}")
    save_posterior_estimates(history, output_dir)  # type: ignore
    visualize_history(history, output_dir)  # type: ignore


def save_posterior_estimates(history: object, output_dir: str) -> None:
    """
    Calculate posterior estimates for mean and median from ABC-SMC run results
    and save them to disk.
    """
    df: pd.DataFrame
    df, w = history.get_distribution()  # type: ignore

    estimates: dict[str, dict[str, float]] = {}
    col: str
    for col in df.columns:
        mean: float = float(np.average(df[col], weights=w))  # type: ignore
        sorted_idx = np.argsort(df[col])  # type: ignore
        cum_weights = np.cumsum(w[sorted_idx])  # type: ignore
        median: float = float(
            df[col].iloc[  # type: ignore
                sorted_idx[np.searchsorted(cum_weights, 0.5)]  # type: ignore
            ]
        )
        estimates[col] = {"mean": mean, "median": median}

    est_df = pd.DataFrame(estimates).T
    est_df.to_csv(os.path.join(output_dir, "abc_estimates.csv"))


def visualize_history(history: object, output_dir: str) -> None:
    """
    Visualize ABC-SMC run results.
    """
    df: pd.DataFrame
    df, w = history.get_distribution()  # type: ignore


def create_argparser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser.

    Returns:
        Configured ArgumentParser instance
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="learn_nuclease_params",
        description="Learn nuclease parameters from input data using ABC-SMC.",
        epilog="{}. Licensed under GPLv3.".format(version_string))

    parser.add_argument(
        "-d", "--motif-data", type=str, dest="motif_data", required=True,
        help="CSV file from pyfrag.py stats command containing 5' end motif "
        "data (e.g., sample_k4_end_motifs.csv)")

    parser.add_argument(
        "-s", "--pop-size", type=int, dest="pop_size", required=False,
        default=10, help="Per-generation population size to use")

    parser.add_argument(
        "-f", "--fasta", type=str, dest="fasta_path", required=True,
        help="Path to reference FASTA file (must be indexed)")

    parser.add_argument(
        "-o", "--output-dir", type=str, dest="output_dir", default=".",
        help="Output directory for results (default: current directory)")

    parser.add_argument(
        "-v", "--verbose", action="store_true", dest="is_verbose",
        default=False, help="Enable verbose logging")

    parser.add_argument(
        "--version", action="version", version=version_string,
        help="Show version information and exit")

    return parser


def main() -> None:
    """
    Entry point for the parameter fitting.
    """
    logger: logging.Logger = logging.getLogger(LOGGER_NAME)
    parser: argparse.ArgumentParser = create_argparser()
    args: argparse.Namespace = parser.parse_args()

    level: int = logging.INFO
    is_verbose: Final[bool] = args.is_verbose
    if is_verbose:
        level = logging.DEBUG
    logger.setLevel(level)

    motif_data_path: str = args.motif_data
    fasta_path: str = args.fasta_path
    pop_size: int = args.pop_size
    output_dir: str = args.output_dir

    if not os.path.exists(motif_data_path):
        fail(f"Motif data file does not exist: {motif_data_path}", logger)

    if not os.path.exists(fasta_path):
        fail(f"FASTA file does not exist: {fasta_path}", logger)

    fasta_index = fasta_path + ".fai"
    if not os.path.exists(fasta_index):
        fail(f"FASTA index not found: {fasta_index}. "
             f"Please run 'samtools faidx {fasta_path}'", logger)

    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            logger.info(f"Created output directory: {output_dir}")
        except OSError as e:
            fail(f"Cannot create output directory: {e}", logger)

    logger.info("Starting nuclease parameter estimation using ABC-SMC")
    logger.info(f"Motif data: {motif_data_path}")
    logger.info(f"Reference FASTA: {fasta_path}")
    logger.info(f"Output directory: {output_dir}")

    try:
        optimize(motif_data_path, fasta_path, output_dir, pop_size, logger)
        logger.info("Parameter estimation completed successfully")
    except Exception as e:
        fail(f"Parameter estimation failed: {e}", logger)


if __name__ == "__main__":
    main()
