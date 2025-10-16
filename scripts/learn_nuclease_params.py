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
import signal
import sys
import pyabc  # type: ignore
import pyfraglib
import warnings

import numpy as np
import numpy.typing as npt
import pandas as pd

from typing import Final, NoReturn
from functools import partial
from pyfraglib import FragmentSimulator, NucleaseProfile, FragmentList
from pyfraglib.simulator.fragment_simulator import SequenceContextGenerator
from pyfraglib.core import detect_cpus
from pyfraglib.stats import end_motifs_barplot

version_string: Final[str] = \
    "learn_nuclease_params v{} (running on Python v{})".format(
        pyfraglib.__version__, sys.version.split(" ")[0]
    )
LOGGER_NAME: Final[str] = "learn_nuclease_params"
EPS: Final[float] = 1e-4
KERNEL_SCALING: Final[float] = 2.0


def bounded_normal_prior(  # type: ignore
    mean: float, std: float, low: float = 0.0, high: float = 2.5
) -> pyabc.RV:
    a, b = (low - mean) / std, (high - mean) / std
    return pyabc.RV(  # type: ignore
        "truncnorm", a=a, b=b, loc=mean, scale=std
    )


class RandomWalkTransition(pyabc.transition.Transition):  # type: ignore
    def __init__(self, step_size: float = 0.01) -> None:
        super().__init__()
        self.step_size: float = step_size

    def rvs(  # type: ignore
        self, source_parameter, source_pdf_val, target_pdf_val, t, **kwargs
    ) -> dict[str, object]:
        """Simple random walk that doesn't depend on other particles."""
        perturbed: dict[str, object] = {}
        for key, val in source_parameter.items():  # type: ignore
            noise = np.random.normal(0, val * self.step_size)  # type: ignore
            perturbed[key] = max(val + noise, EPS)  # type: ignore
            return perturbed

    def pdf(  # type: ignore
        self, source_parameter, target_parameter, source_pdf_val,
        target_pdf_val, t
    ) -> float:
        """Return constant probability."""
        return 1.0


def fail(
    msg: str, logger: logging.Logger | None
) -> NoReturn:
    """
    Log a fatal error message and exit the program.

    Kills all processes of the process group, too.

    Args:
        msg: Error message to log
        logger: Logger instance to use for logging

    Raises:
        SystemExit: Always exits with code 1
    """
    if logger:
        logger.fatal(msg)
    else:
        print(msg)
    os.killpg(os.getpgid(os.getpid()), signal.SIGTERM)
    sys.exit(1)


def kill_on_warning(
    message: str, _category: object, _filename: str,
    _lineno: int, file: object | None = None,
    line: object | None = None,
    logger: logging.Logger | None = None
) -> NoReturn:
    fail(f"failed with warning: {message}", logger)


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
        motif: (count / total_count) * 1000
        for motif, count in motif_counts.items()
    }
    return motif_freqs


def simulate_fragments(
    params: dict[str, float], fasta_path: str, num_fragments: int
) -> FragmentList:
    """
    Simulate fragments given nuclease parameters.

    Args:
        params: Dict containing nuclease activity and preference parameters
        fasta_path: path to reference fasta for sequence generation
        num_fragments: number of fragments to generate

    Returns:
        FragmentList
    """
    dnase1_cg_pref: float = max(params["dnase1_cg_pref"], EPS)
    dnase1_at_pref: float = max(params["dnase1_at_pref"], EPS)
    dnase1_aa_pref: float = max(params["dnase1_aa_pref"], EPS)
    dnase1_motif_preference: dict[str, float] = {
        "CG": dnase1_cg_pref,
        "GC": dnase1_cg_pref,
        "AT": dnase1_at_pref,
        "TA": dnase1_at_pref,
        "AA": dnase1_aa_pref,
        "TT": dnase1_aa_pref
    }

    dnase1l3_cc_pref: float = max(params["dnase1l3_cc_pref"], EPS)
    dnase1l3_c_pref: float = max(params["dnase1l3_c_pref"], EPS)
    dnase1l3_cg_pref: float = max(params["dnase1l3_cg_pref"], EPS)
    dnase1l3_ct_pref: float = max(params["dnase1l3_ct_pref"], EPS)
    dnase1l3_gg_pref: float = max(params["dnase1l3_gg_pref"], EPS)
    dnase1l3_at_pref: float = max(params["dnase1l3_at_pref"], EPS)
    dnase1l3_a_pref: float = max(params["dnase1l3_a_pref"], EPS)
    dnase1l3_motif_preference: dict[str, float] = {
        "CC": dnase1l3_cc_pref,
        "C": dnase1l3_c_pref,
        "CG": dnase1l3_cg_pref,
        "GC": dnase1l3_cg_pref,
        "CT": dnase1l3_ct_pref,
        "TC": dnase1l3_ct_pref,
        "GG": dnase1l3_gg_pref,
        "AT": dnase1l3_at_pref,
        "TA": dnase1l3_at_pref,
        "A": dnase1l3_a_pref,
        "T": dnase1l3_a_pref
    }

    dffb_a_pref: float = max(params["dffb_a_pref"], EPS)
    dffb_c_pref: float = max(params["dffb_c_pref"], EPS)
    dffb_motif_preference: dict[str, float] = {
        "A": dffb_a_pref,
        "T": dffb_a_pref,
        "C": dffb_c_pref,
        "G": dffb_c_pref
    }

    dnase1_activity: float = max(params["dnase1_activity"], EPS)
    dnase1l3_activity: float = max(params["dnase1l3_activity"], EPS)
    dffb_activity: float = max(params["dffb_activity"], EPS)
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
            num_fragments=num_fragments, tissue_type="healthy",
            nuclease_profile=nuclease_profile
        )
        return fragments

    finally:
        seq_gen.close()


def simulate_motifs(
    params: dict[str, float], fasta_path: str, num_fragments: int
) -> dict[str, float]:
    """
    Simulate fragment 5' end motifs given nuclease parameters.

    Args:
        params: Dict containing nuclease activity and preference parameters
        fasta_path: path to reference fasta for sequence generation
        num_fragments: number of fragments to generate

    Returns:
        Dictionary of end motif frequencies
    """
    fragments: FragmentList = simulate_fragments(
        params, fasta_path, num_fragments
    )
    motifs_5p, _, num_frags = fragments.count_endmotifs(kmer_len=4)
    if num_frags == 0:
        raise ValueError("No 5' end motifs simulated")

    total_count: int = sum(motifs_5p.values())
    motif_freqs: dict[str, float] = {
        motif: (count / total_count) * 1000
        for motif, count in motifs_5p.items()
    }
    return motif_freqs


def calculate_distance(
    simulated: dict[str, float], observed: dict[str, float],
    dist: str = "L1"
) -> float:
    """
    Calculate distance between simulated and observed motif frequencies.

    Note:
        The only difference to L2 is that we skip taking the square root.

    Args:
        simulated: Simulated motif frequencies
        observed: Observed motif frequencies from data
        dist: either "L1" or "L2" (only L2-like without sqrt)

    Returns:
        The requested distance metric
    """
    all_motifs = set(simulated.keys()) | set(observed.keys())
    distance: float = 0.0
    for motif in all_motifs:
        sim_freq = simulated.get(motif, 0.0)
        obs_freq = observed.get(motif, 0.0)
        if dist == "L1":
            distance += abs(sim_freq - obs_freq)
        elif dist == "L2":
            distance += (sim_freq - obs_freq) ** 2
        else:
            raise ValueError(f"unknown distance {dist} requested")
    return max(distance, EPS)


def abc_model(
    params: dict[str, float], fasta_path: str, num_fragments: int
) -> dict[str, float]:
    """
    ABC model function that simulates motifs given parameters.
    """
    return simulate_motifs(params, fasta_path, num_fragments)


def abc_distance(
    simulated: dict[str, float], observed: dict[str, float]
) -> float:
    """
    ABC distance function that compares simulated to observed data.
    """
    return calculate_distance(simulated, observed)


def optimize(
    motif_data_path: str, fasta_path: str, output_dir: str,
    num_fragments: int, max_pops: int, pop_size: int,
    logger: logging.Logger
) -> None:
    """
    Run ABC-SMC parameter estimation.

    Args:
        motif_data_path: Path to CSV file with observed motif data
        fasta_path: Path to reference FASTA file
        output_dir: Directory for output files
        num_fragments: Number of fragments per simulation run
        max_pops: Maximum number of populations / iterations
        pop_size: Population size per iteration
        logger: Logger instance
    """
    observed_motifs = load_observed_motifs(motif_data_path, logger)
    logger.info(f"Loaded {len(observed_motifs)} motifs from observed data")

    prior = pyabc.Distribution(  # type: ignore
        dnase1_activity=(
            pyabc.RV("uniform", loc=0.0, scale=2.0)  # type: ignore
        ),
        dnase1l3_activity=(
            pyabc.RV("uniform", loc=0.0, scale=2.0)  # type: ignore
        ),
        dffb_activity=(
            pyabc.RV("uniform", loc=0.0, scale=2.0)  # type: ignore
        ),

        dnase1_cg_pref=(
            bounded_normal_prior(0.8, 0.4)  # type: ignore
        ),
        dnase1_at_pref=(
            bounded_normal_prior(1.2, 0.4)  # type: ignore
        ),
        dnase1_aa_pref=(
            bounded_normal_prior(1.1, 0.4)  # type: ignore
        ),

        dnase1l3_cc_pref=(
            bounded_normal_prior(6.0, 2.0, high=12.0)  # type: ignore
        ),
        dnase1l3_c_pref=(
            bounded_normal_prior(3.0, 1.5, high=10.0)  # type: ignore
        ),
        dnase1l3_cg_pref=(
            bounded_normal_prior(2.5, 1.4, high=8.0)  # type: ignore
        ),
        dnase1l3_ct_pref=(
            bounded_normal_prior(2.0, 1.3, high=6.0)  # type: ignore
        ),
        dnase1l3_gg_pref=(
            bounded_normal_prior(1.8, 1.2, high=5.0)  # type: ignore
        ),
        dnase1l3_at_pref=(
            bounded_normal_prior(0.6, 0.2)  # type: ignore
        ),
        dnase1l3_a_pref=(
            bounded_normal_prior(0.7, 0.2)  # type: ignore
        ),

        dffb_a_pref=(
            bounded_normal_prior(1.1, 0.5)  # type: ignore
        ),
        dffb_c_pref=(
            bounded_normal_prior(0.9, 0.5)  # type: ignore
        )
    )

    n_cores = detect_cpus()
    logger.info(f"Using {n_cores} cores for ABC-SMC sampling")
    sampler = pyabc.sampler.MulticoreEvalParallelSampler(  # type: ignore
        n_procs=n_cores
    )

    transition = pyabc.MultivariateNormalTransition(  # type: ignore
        scaling=KERNEL_SCALING
    )
    # @NOTE(ds): A local transition kernel is a valid alternative:
    #
    # transition = pyabc.LocalTransition(  # type: ignore
    #     scaling=0.01, k_fraction=0.8
    # )

    abc = pyabc.ABCSMC(  # type: ignore
        models=partial(
            abc_model, fasta_path=fasta_path, num_fragments=num_fragments
        ),
        parameter_priors=prior,  # type: ignore
        distance_function=abc_distance,
        population_size=pop_size,
        transitions=transition,  # type: ignore
        sampler=sampler,  # type: ignore
        eps=pyabc.QuantileEpsilon(alpha=0.9)  # type: ignore
    )

    db_path: str = f"sqlite:///{output_dir}/abc_history.db"
    history: pyabc.History = abc.new(  # type: ignore
        db_path, observed_motifs
    )

    logger.info(f"Running ABC-SMC inference with DB at {db_path}")
    history = abc.run(  # type: ignore
        minimum_epsilon=0.01,
        max_nr_populations=max_pops,
        min_acceptance_rate=0.01
    )

    logger.info(
        f"Completed generations: {history.max_t + 1}"  # type: ignore
    )


def analyze_run_from_history(
    db_file: str, output_dir: str, fasta_path: str, logger: logging.Logger
) -> None:
    """
    Analyze an existing run.
    """
    logger.info(f"loading history from file `{db_file}`")
    db_str: str = f"sqlite:///{db_file}"
    history: pyabc.History = pyabc.History(db_str)  # type: ignore

    df: pd.DataFrame
    w: npt.NDArray[np.float64]
    df, w = history.get_distribution()  # type: ignore

    logger.info(f"saving posterior estimates to `{output_dir}`")
    estimates_df = save_posterior_estimates(df, w, output_dir)

    logger.info("performing simulation based on posterior estimates")
    nuc_params: dict[str, float] = dict(estimates_df["median"])  # type: ignore
    fragments: FragmentList = simulate_fragments(
        nuc_params, fasta_path, 100_000
    )
    end_motifs_barplot(fragments, output_dir, "simulation_after_fit", 4)


def save_posterior_estimates(
    df: pd.DataFrame, w: npt.NDArray[np.float64], output_dir: str
) -> pd.DataFrame:
    """
    Calculate posterior estimates for mean and median from ABC-SMC run results
    and save them to disk.
    """
    estimates: dict[str, dict[str, float]] = {}
    param: str
    for param in df.columns:
        mean: float = float(np.average(df[param], weights=w))  # type: ignore
        std: float = float(
            np.sqrt(  # type: ignore
                np.average((df[param] - mean)**2, weights=w)  # type: ignore
            )
        )
        sorted_indices = np.argsort(df[param])  # type: ignore
        sorted_values = df[param].values[sorted_indices]  # type: ignore
        sorted_weights = w[sorted_indices]  # type: ignore
        cumsum = np.cumsum(sorted_weights)  # type: ignore

        estimates[param] = {
            "mean": mean,
            "std": std,
            "median":
                sorted_values[np.searchsorted(cumsum, 0.5)],  # type: ignore
            "min": float(df[param].min()),  # type: ignore
            "max": float(df[param].max()),  # type: ignore
        }

    est_df = pd.DataFrame(estimates).T
    est_df.to_csv(os.path.join(output_dir, "abc_estimates.csv"))

    return est_df


def visualize_history(
    df: pd.DataFrame, w: npt.NDArray[np.float64], output_dir: str
) -> None:
    """
    Visualize ABC-SMC run results.
    """
    pass


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
        "-d", "--motif-data", type=str, dest="motif_data", required=False,
        default=None, help="CSV file from pyfrag.py stats command containing "
        "5' end motif data (e.g., sample_k4_end_motifs.csv)")

    parser.add_argument(
        "-n", "--num-fragments", type=int, dest="num_fragments",
        required=False, default=10_000,
        help="Number of fragments per simulation run")

    parser.add_argument(
        "-s", "--pop-size", type=int, dest="pop_size", required=False,
        default=10, help="Per-generation population size to use")

    parser.add_argument(
        "-p", "--max-pops", type=int, dest="max_pops", required=False,
        default=15, help="Maximum number of populations")

    parser.add_argument(
        "-f", "--fasta", type=str, dest="fasta_path", required=True,
        help="Path to reference FASTA file (must be indexed)")

    parser.add_argument(
        "-o", "--output-dir", type=str, dest="output_dir", default=".",
        help="Output directory for results (default: current directory)")

    parser.add_argument(
        "--analyze-run", action="store_true", dest="analyze_run",
        default=False, help="Instead of learning parameters, analyze a "
        "past run (requires --db-file to be set)")

    parser.add_argument(
        "--db-file", type=str, dest="db_file", required=False, default=None,
        help="Indicates the database file (only used when called "
        "with --analyze-run")

    parser.add_argument(
        "-v", "--verbose", action="store_true", dest="is_verbose",
        default=False, help="Enable verbose logging")

    parser.add_argument(
        "--version", action="version", version=version_string,
        help="Show version information and exit")

    return parser


def main() -> None:
    """
    Entry point for the parameter fitting routine.
    """
    pyabc_logger: logging.Logger = logging.getLogger("ABC")
    pyabc_logger.handlers.clear()
    pyabc_logger.propagate = True

    logger: logging.Logger = logging.getLogger(LOGGER_NAME)
    parser: argparse.ArgumentParser = create_argparser()
    args: argparse.Namespace = parser.parse_args()

    level: int = logging.INFO
    is_verbose: Final[bool] = args.is_verbose
    if is_verbose:
        level = logging.DEBUG
    logger.setLevel(level)

    warnings.showwarning = partial(
        kill_on_warning, logger=logger
    )

    fasta_path: str = args.fasta_path
    if not fasta_path or not os.path.exists(fasta_path):
        fail(f"FASTA file does not exist: {fasta_path}", logger)

    fasta_index: str = fasta_path + ".fai"
    if not os.path.exists(fasta_index):
        fail(f"FASTA index not found: {fasta_index}. "
             f"Please run 'samtools faidx {fasta_path}'", logger)

    output_dir: str = args.output_dir
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            logger.info(f"Created output directory: {output_dir}")
        except OSError as e:
            fail(f"Cannot create output directory: {e}", logger)

    analyze_run: bool = args.analyze_run
    if analyze_run:
        db_file: str | None = args.db_file
        if not db_file or not os.path.exists(db_file):
            fail(f"DB file does not exist: {db_file}", logger)
        analyze_run_from_history(
            db_file, output_dir, fasta_path, logger
        )
        return

    motif_data_path: str | None = args.motif_data
    num_fragments: int = args.num_fragments
    max_pops: int = args.max_pops
    pop_size: int = args.pop_size

    if not motif_data_path or not os.path.exists(motif_data_path):
        fail(f"Motif data file does not exist: {motif_data_path}", logger)

    logger.info("Starting nuclease parameter estimation using ABC-SMC")
    logger.info(f"Motif data: {motif_data_path}")
    logger.info(f"Reference FASTA: {fasta_path}")
    logger.info(f"Output directory: {output_dir}")

    try:
        optimize(
            motif_data_path, fasta_path, output_dir,
            num_fragments, max_pops, pop_size, logger
        )
        logger.info("Parameter estimation completed successfully")
    except Exception as e:
        fail(f"Parameter estimation failed: {e}", logger)


if __name__ == "__main__":
    main()
