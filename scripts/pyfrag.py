#!/usr/bin/env python3
"""
PyFragLib Command Line Interface
================================

This module provides a command-line interface for the PyFragLib software suite,
which calculates fragmentomics features from cfDNA and performs downstream
analyses.

The CLI supports multiple subcommands for different analysis workflows:
- extract: Extract fragments from BAM files
- stats: Generate descriptive statistics from fragment files
- lengths: Analyze fragment length distributions with Gaussian mixture models
- scores: Calculate fragmentomics scores (WPS, motif diversity)
- simulate: Generate synthetic cfDNA fragments based on biological parameters

Each subcommand has specific options and requirements. Use
`pyfrag.py <subcommand> -h`
for detailed help on individual subcommands.

Examples:
    Extract fragments from a BAM file:
        pyfrag.py extract --bam-file sample.bam --out-dir output/

    Generate statistics from fragment files:
        pyfrag.py stats --frag-file sample.frag --out-dir results/

    Simulate cfDNA fragments:
        pyfrag.py simulate --config simulation.json --out-dir sim_output/

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
import signal
import sys
import pysam
import pyfraglib

import pandas as pd

from functools import partial
from multiprocessing import Pool
from typing import Final, NoReturn, Optional
from pyfraglib import Fragment, FragmentList
from pyfraglib.core import CodeUnreachableError, parse_bed_file, detect_cpus
from pyfraglib.fragfile import FragFile
from pyfraglib.lengths import fragment_length_plot, fragment_length_gmm
from pyfraglib.stats import fragments_per_chromosome_barplot, \
                            end_motifs_barplot, log_stats, \
                            export_length_distribution_csv, \
                            export_end_motifs_csv
from pyfraglib.scores import motif_diversity, windowed_protection_score, \
                             score_line_plot
from pyfraglib.simulator import FragmentSimulator, NucleaseProfile, \
                                TissueMixtureSimulator, \
                                SequenceContextGenerator

version_string: Final[str] = "pyfraglib v{} (running on Python v{})" \
    .format(pyfraglib.__version__, sys.version.split(" ")[0])


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


def signal_handler(sig: int, frame: object) -> NoReturn:
    """
    Handle interrupt signals (e.g., Ctrl+C) gracefully.

    Args:
        sig: Signal number
        frame: Current stack frame

    Raises:
        SystemExit: Always exits with code 1
    """
    logger: logging.Logger = logging.getLogger("pyfrag")
    fail("an error occurred", logger)


def create_argparser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser for the CLI.

    Sets up the main argument parser with global options and subparsers for
    each subcommand (extract, stats, lengths, scores, simulate, version).

    Returns:
        Configured ArgumentParser instance with all subcommands

    Note:
        We always require our users to indicate an output directory. What
        type of output we produce is determined by a subcommand (implemented
        via argparse's subparsers). Subcommands have a bunch of additional, but
        individual options which determine their behavior only.
    """
    parent_parser: argparse.ArgumentParser = \
        argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        "-o", "--out-dir", type=str, dest="out_dir", default="pyfrag_out",
        help="Directory that output is saved to. The type of output depends "
        "on the subcommand. Dir is created if it does not exist.")
    parent_parser.add_argument(
        "-v", "--verbose", action="store_true", dest="is_verbose",
        default=False, help="Set log level so that debugging info is printed.")

    argparser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="pyfrag", description="Use a subset of `pyfraglib's "
        "capabilities from the command line.",
        epilog="{}. Licensed under GPLv3. See repository at `{}' for more "
        "info.".format(version_string, pyfraglib.__repo_url__))
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser] = \
        argparser.add_subparsers(
            required=True, help="Select from a number of utilities. "
            "Example: ``$ pyfraglib.py extract --bam-file=<FILE>''",
            dest="subcommand")

    subparsers.add_parser("version", help="Show version info.",
        parents=[parent_parser]
    )

    argparser_extract: argparse.ArgumentParser = subparsers.add_parser(
        "extract", help="Extract fragment information from BAM file(s). "
        "Results are written to disk in the `.frag' file format.",
        parents=[parent_parser]
    )
    argparser_extract.add_argument(
        "-d", "--bam-dir", type=str, dest="bam_dir", required=False,
        help="Input BAM files to be analyzed in batch mode. Must not be "
        "combined with `--bam-file'. Might consume a lot of memory.")
    argparser_extract.add_argument(
        "-f", "--bam-file", type=str, dest="bam_file", required=False,
        help="Input BAM file to be analyzed. Must not be combined with "
        "`--bam-dir'.")
    argparser_extract.add_argument(
        "--with-vcf", action="store_true", dest="has_vcf_file",
        default=False, help="If flag is set, there must be a VCF file for "
        "every BAM (e.g. `patient1.bam' must have `patient1.vcf' in the same "
        "directory). Variants must be SNVs. All fragments for which one of "
        "the reads carries the ALT allele are then annotated as mutated.")
    argparser_extract.add_argument(
        "--nanopore", action="store_true", dest="is_nanopore",
        default=False, help="If set, input BAM file(s) will be treated as "
        "single-end long read sequencing (e.g. ONT sequencing). If not set "
        "(default), BAMs are expected to be paired-end.")

    argparser_stats: argparse.ArgumentParser = subparsers.add_parser(
        "stats", help="Given (a) `.frag' file(s), extract descriptive "
        "statistics and save them to disk.", parents=[parent_parser])
    argparser_stats.add_argument(
        "-d", "--frag-dir", type=str, dest="frag_dir", required=False,
        help="Input FRAG files to be analyzed in batch mode. Must not be "
        "combined with `--frag-file'. Requires sufficient amounts of "
        "memory.")
    argparser_stats.add_argument(
        "-f", "--frag-file", type=str, dest="frag_file", required=False,
        help="Input FRAG file to be analyzed. Must not be combined with "
        "`--frag-dir'.")
    argparser_stats.add_argument(
        "--kmer-length", type=int, dest="kmer_len", default=4,
        help="K-mer length for end motif analyses."
        )

    argparser_lengths: argparse.ArgumentParser = subparsers.add_parser(
        "lengths", help="Given (a) `.frag' file(s), read fragment length "
        "data and create a number of plots that are saved to disk.",
        parents=[parent_parser])
    argparser_lengths.add_argument(
        "-d", "--frag-dir", type=str, dest="frag_dir", required=False,
        help="Input FRAG files to be analyzed in batch mode. Must not be "
        "combined with `--frag-file'. Requires sufficient amounts of "
        "memory.")
    argparser_lengths.add_argument(
        "-f", "--frag-file", type=str, dest="frag_file", required=False,
        help="Input FRAG file to be analyzed. Must not be combined with "
        "`--frag-dir'.")
    argparser_lengths.add_argument(
        "-c", "--config-file", type=str, dest="config_file", required=True,
        help="JSON config file for Gaussian mixture model. See "
        "`config/gmm.json' for an example.")

    argparser_scores: argparse.ArgumentParser = subparsers.add_parser(
        "scores", help="Given (a) `.frag' file(s), calculate a variety of "
        "fragmentomics scores and save the results to disk.",
        parents=[parent_parser])
    argparser_scores.add_argument(
        "-d", "--frag-dir", type=str, dest="frag_dir", required=False,
        help="Input FRAG files to be analyzed in batch mode. Must not be "
        "combined with `--frag-file'. Requires sufficient amounts of "
        "memory.")
    argparser_scores.add_argument(
        "-f", "--frag-file", type=str, dest="frag_file", required=False,
        help="Input FRAG file to be analyzed. Must not be combined with "
        "`--frag-dir'.")
    argparser_scores.add_argument(
        "-b", "--bed-file", type=str, dest="bed_file", required=True,
        help="A correctly formatted BED file (see docs) that provides "
        "genomic coordinates for score calculation.")
    argparser_scores.add_argument(
        "--hg38", action="store_true", dest="is_hg38",
        default=False, help="If set, hg38 will be assumed for the analyses. "
        "If unset (default), hg19 is assumed.")

    argparser_simulate: argparse.ArgumentParser = subparsers.add_parser(
        "simulate", help="Simulate synthetic cfDNA fragments based on "
        "biological parameters. Results are written to a `.frag' file.",
        parents=[parent_parser])
    argparser_simulate.add_argument(
        "-c", "--config", type=str, dest="config_file", required=True,
        help="JSON configuration file containing simulation parameters. "
        "See documentation for configuration format.")

    return argparser


def extract(out_dir: str, args: argparse.Namespace) -> None:
    """
    Extract fragment information from BAM files.

    Processes BAM files to extract cfDNA fragments and saves them in the
    .frag format. Supports both single-file and batch processing modes.
    Can handle paired-end and single-end (Nanopore) sequencing data.

    Args:
        out_dir: Output directory for .frag files
        args: Command line arguments containing BAM file paths and options

    Raises:
        SystemExit: If invalid arguments or file not found

    Note:
        Batch processing of multiple BAM files is no longer recommended. The
        preferred way to handle large datasets is using the Nextflow pipeline
        included with PyFragLib.

        At this point, we can assume that `out_dir' exists. Loading all
        `FragmentList's into a collection before saving the to disk is super
        expensive in terms of memory. We did that in previous versions of
        `pyfraglib' but changed things around to be able to handle large
        directories of even larger BAM files.
    """
    bam_file: Final[str] = args.bam_file
    bam_dir: Final[str] = args.bam_dir
    is_nanopore: Final[bool] = args.is_nanopore

    if bam_file and bam_dir:
        fail("--bam-dir and --bam-file are incompatible options", logger)
    elif not bam_file and not bam_dir:
        fail("either one of --bam-dir or --bam-file is required", logger)

    bam_files: list[str] = []
    if bam_dir:
        if not os.path.isdir(bam_dir):
            fail("directory `{}' does not exist".format(bam_dir), logger)
        logger.info("looking for BAM files in `{}'".format(bam_dir))
        bam_files = search_dir(bam_dir, [".BAM", ".bam"])
    elif bam_file:
        if not os.path.isfile(bam_file):
            fail("file `{}' does not exist".format(bam_file), logger)
        bam_files.append(bam_file)
    else:
        raise CodeUnreachableError("internal error (no bam file or dir)")

    if len(bam_files) == 0:
        fail("no bam files found; file extension `.bam' is required", logger)

    has_vcf_file: Final[bool] = args.has_vcf_file
    vcf_files: Optional[list[str]] = None
    if has_vcf_file:
        vcf_files = []
        for bam_file in bam_files:
            basename, _ = os.path.splitext(bam_file)
            vcf_name = "{}.vcf".format(basename)
            vcf_files.append(vcf_name)
            if not os.path.isfile(vcf_name):
                fail("VCF file `{}' not found".format(vcf_name), logger)
            logger.info("found VCF file `{}'".format(vcf_name))
    else:
        logger.info("not expecting any VCF files")

    # NOTE(ds): Instead of accumulating all results in memory and writing them
    # out at once, we extract and write BAM/FRAG files independently.
    Fragment.bams_to_frags(bam_files, vcf_files, out_dir, is_nanopore)


def stats(out_dir: str, args: argparse.Namespace) -> None:
    """
    Generate descriptive statistics from fragment files.

    Analyzes .frag files to produce various statistics and visualizations:
    - Fragment count per chromosome
    - End motif distributions
    - Summary statistics logged to console

    Args:
        out_dir: Output directory for statistics files and plots
        args: Command line arguments containing fragment file paths

    Raises:
        SystemExit: If invalid arguments or file not found
    """
    frag_file: Final[str] = args.frag_file
    frag_dir: Final[str] = args.frag_dir
    kmer_len: Final[int] = args.kmer_len

    if frag_file and frag_dir:
        fail("--frag-dir and --frag-file are incompatible options", logger)
    elif not frag_file and not frag_dir:
        fail("either one of --frag-dir or --frag-file is required", logger)

    frag_files: list[str] = get_frag_files(frag_dir, frag_file)
    if len(frag_files) == 0:
        fail("no bam files found; file extension `.frag' is required", logger)

    for path in frag_files:
        logger.info("loading `{}'".format(path))

        fragment_file: FragFile = FragFile(path)
        fragments: FragmentList = fragment_file.get_fragment_list()
        fragment_file.close()

        name: str = filename_only(path)
        log_stats(fragments, logger, out_dir, name)
        fragments_per_chromosome_barplot(fragments, out_dir, name)
        export_length_distribution_csv(fragments, out_dir, name)

        end_motifs_barplot(fragments, out_dir, name, kmer_len=kmer_len)
        export_end_motifs_csv(fragments, out_dir, name, kmer_len=kmer_len)


def lengths(out_dir: str, args: argparse.Namespace) -> None:
    """
    Analyze fragment length distributions.

    Processes fragment files to analyze length distributions using:
    - Fragment length histograms and plots
    - Gaussian mixture model fitting based on configuration

    Args:
        out_dir: Output directory for length analysis results
        args: Command line arguments containing fragment files and config

    Raises:
        SystemExit: If invalid arguments or missing config file
    """
    frag_file: Final[str] = args.frag_file
    frag_dir: Final[str] = args.frag_dir
    config_file: Final[str] = args.config_file

    if frag_file and frag_dir:
        fail("--frag-dir and --frag-file are incompatible options", logger)
    elif not frag_file and not frag_dir:
        fail("either one of --frag-dir or --frag-file is required", logger)

    if not os.path.exists(config_file):
        fail("config file `{}' does not exist".format(config_file), logger)

    frag_files: list[str] = get_frag_files(frag_dir, frag_file)
    if len(frag_files) == 0:
        fail("no bam files found; file extension `.frag' is required", logger)

    for path in frag_files:
        logger.info("loading `{}'".format(path))

        fragment_file: FragFile = FragFile(path)
        fragments: FragmentList = fragment_file.get_fragment_list()
        fragment_file.close()

        name: str = filename_only(path)
        fragment_length_plot(fragments, out_dir, name)
        fragment_length_gmm(fragments, config_file, out_dir, name)


def scores(out_dir: str, args: argparse.Namespace) -> None:
    """
    Calculate fragmentomics scores from fragment files.

    Computes various fragmentomics metrics including:
    - Windowed Protection Score (WPS)
    - Motif diversity scores (Shannon entropy, Simpson index)
    - Score visualization plots

    Args:
        out_dir: Output directory for score results and plots
        args: Command line arguments containing fragment files and BED regions

    Raises:
        SystemExit: If invalid arguments or missing required files

    Note:
        We interleave the calculation of our scores. It's ugly, but otherwise
        we have to store all fragment files as a collection which requires
        loads of memory.
    """
    frag_file: Final[str] = args.frag_file
    frag_dir: Final[str] = args.frag_dir
    bed_file: Final[str] = args.bed_file

    is_hg38: Final[bool] = args.is_hg38
    genome: Final[str] = "hg38" if is_hg38 else "hg19"

    if frag_file and frag_dir:
        fail("--frag-dir and --frag-file are incompatible options", logger)
    elif not frag_file and not frag_dir:
        fail("either one of --frag-dir or --frag-file is required", logger)

    frag_files: list[str] = get_frag_files(frag_dir, frag_file)
    if len(frag_files) == 0:
        fail("no bam files found; file extension `.frag' is required", logger)

    regions: pysam.TabixFile = pysam.TabixFile(bed_file)
    col_names: list[str] = ["sample_name", "shannon_3p", "simpson_3p",
                            "shannon_5p", "simpson_5p"]
    glbl_motif_diversity_scores: pd.DataFrame = pd.DataFrame(columns=col_names)

    for it, path in enumerate(frag_files):
        logger.info("loading `{}'".format(path))

        fragment_file: FragFile = FragFile(path)
        fragments: FragmentList = fragment_file.get_fragment_list()
        fragment_file.close()

        name: str = filename_only(path)
        shannon_5p, shannon_3p = motif_diversity(fragments, name, "shannon")
        simpson_5p, simpson_3p = motif_diversity(fragments, name, "simpson")

        new_row: list[str | float] = [
            name, shannon_3p, simpson_3p, shannon_5p, simpson_5p
        ]
        glbl_motif_diversity_scores.loc[it] = new_row

        logger.info("calculating windowed protection score")
        wps_df: pd.DataFrame = windowed_protection_score(fragments, regions,
                                                         genome=genome)
        wps_outpath: str = os.path.join(out_dir, "wps_{}.csv".format(name))
        logger.info("saving windowed protection score to `{}'".format(
            wps_outpath))
        wps_df.to_csv(wps_outpath, index=False)

        score_line_plot(wps_df, name, out_dir, genome=genome)

    glbl_mds_outpath: str = os.path.join(
        out_dir, "global_motif_diversity_scores.csv"
    )
    logger.info("saving global motif diversity scores for all files "
                "to `{}'".format(glbl_mds_outpath))
    glbl_motif_diversity_scores.to_csv(glbl_mds_outpath, index=False)


def simulate(out_dir: str, args: argparse.Namespace) -> None:
    """
    Simulate synthetic cfDNA fragments based on configuration file.

    Supports multiple simulation modes:
    - Basic simulation: single tissue type with specified parameters
    - Tissue mixture: multiple tissues with specified fractions
    - Cancer progression: tumor fraction changes over time

    Args:
        out_dir: Output directory for .frag file(s)
        args: Command line arguments containing config_file
    """
    config_file: Final[str] = args.config_file
    if not os.path.exists(config_file):
        fail("config file `{}' does not exist".format(config_file), logger)

    logger.info("loading configuration from `{}'".format(config_file))
    with open(config_file, "r") as f:
        config: dict[str, object] = json.load(f)

    try:
        output_name: str = config["output_name"]  # type: ignore
        fasta_path: str = config["fasta_path"]  # type: ignore
        simulation_mode: str = config.get(
            "simulation_mode", "basic"
        )  # type: ignore

        if not os.path.exists(fasta_path):
            fail(
                "FASTA file `{}' does not exist".format(fasta_path), logger
            )

        if simulation_mode == "basic":
            simulate_basic(config, output_name, fasta_path, out_dir)
        elif simulation_mode == "tissue_mixture":
            simulate_tissue_mixture(config, output_name, fasta_path, out_dir)
        elif simulation_mode == "cancer_progression":
            simulate_cancer_progression(
                config, output_name, fasta_path, out_dir
            )
        else:
            fail(
                "unknown simulation mode: {}".format(simulation_mode), logger
            )

    except KeyError as e:
        fail(
            "missing required configuration parameter: {}".format(e), logger
        )


def simulate_region(
    region_data: tuple[str, int, int], fasta_path: str,
    tissue_type: str, nuclease_profile: "NucleaseProfile",
    fragment_size_params: dict[str, float] | None, fragments_per_region: int
) -> "FragmentList":
    """Simulate fragments for a single genomic region."""
    chrom, start, end = region_data
    if fragments_per_region == 0:
        return FragmentList()

    seq_gen = SequenceContextGenerator(fasta_path)
    simulator: FragmentSimulator = FragmentSimulator(seq_gen)
    return simulator.simulate_fragments(
        chrom=chrom, start=start, end=end,
        num_fragments=fragments_per_region,
        tissue_type=tissue_type, nuclease_profile=nuclease_profile,
        fragment_size_params=fragment_size_params
    )


def simulate_basic(
    config: dict[str, object], output_name: str, fasta_path: str, out_dir: str
) -> None:
    """Basic single-tissue simulation."""
    try:
        bed_file: str = config["bed_file"]  # type: ignore
        total_fragments: int = config.get(
            "total_fragments", 10000
        )  # type: ignore
        fragment_params: dict[str, object] = config.get(
            "fragment_params", {}
        )  # type: ignore
        nuclease_params: dict[str, object] = config.get(
            "nuclease_params", {}
        )  # type: ignore
        nuclease_profile: NucleaseProfile = NucleaseProfile(
            dnase1_activity=nuclease_params.get(
                "dnase1_activity", 1.0
            ),  # type: ignore
            dnase1l3_activity=nuclease_params.get(
                "dnase1l3_activity", 1.0
            ),  # type: ignore
            dffb_activity=nuclease_params.get(
                "dffb_activity", 1.0
            )  # type: ignore
        )
        tissue_type: str = config.get(
            "tissue_type", "healthy"
        )  # type: ignore

        if not os.path.exists(bed_file):
            fail("BED file '{}' does not exist".format(bed_file), logger)

        genomic_regions: list[tuple[str, int, int]] = parse_bed_file(bed_file)
    except KeyError as e:
        fail(
            "missing required parameter for basic simulation: {}".format(e),
            logger
        )

    logger.info(
        "running basic simulator with FASTA: {}".format(fasta_path)
    )
    all_fragments: FragmentList = FragmentList()
    fragments_per_region = total_fragments // len(genomic_regions)
    fragment_size_params: dict[str, float] | None = None
    if fragment_params:
        fragment_size_params = {}
        for k, v in fragment_params.items():
            if isinstance(v, (int, float, str)):
                fragment_size_params[k] = float(v)

    logger.info(
        "simulating {} fragments across {} regions in parallel".format(
            total_fragments, len(genomic_regions)
        )
    )

    with Pool(processes=detect_cpus()) as pool:
        partial_simulate = partial(
            simulate_region, fasta_path=fasta_path, tissue_type=tissue_type,
            nuclease_profile=nuclease_profile,
            fragment_size_params=fragment_size_params,
            fragments_per_region=fragments_per_region
        )
        region_results = pool.map(partial_simulate, genomic_regions)

    for fragment_list in region_results:
        for fragment in fragment_list:
            all_fragments.append(fragment)
    logger.info(
        "saving {} simulated fragments to {}.frag".format(
            all_fragments.length(), output_name
        )
    )
    all_fragments.to_frag_file(output_name, out_dir)


def simulate_tissue_mixture(
    config: dict[str, object], output_name: str, fasta_path: str, out_dir: str
) -> None:
    """Tissue mixture simulation."""
    try:
        tissue_types: list[str] = config["tissue_types"]  # type: ignore
        tissue_fractions: list[float] = config[
            "tissue_fractions"
        ]  # type: ignore
        total_fragments: int = config["total_fragments"]  # type: ignore
        bed_file: str = config["bed_file"]  # type: ignore
        add_noise: bool = config.get("add_noise", True)  # type: ignore

        if not os.path.exists(bed_file):
            fail("BED file '{}' does not exist".format(bed_file), logger)

        genomic_regions: list[tuple[str, int, int]] = parse_bed_file(bed_file)

    except KeyError as e:
        fail(
            "missing required parameter for tissue mixture: {}".format(e),
            logger
        )

    logger.info(
        "initializing tissue mixture simulator with FASTA: {}".format(
            fasta_path
        )
    )
    seq_gen: SequenceContextGenerator = SequenceContextGenerator(fasta_path)
    simulator: TissueMixtureSimulator = TissueMixtureSimulator(seq_gen)

    logger.info(
        "simulating tissue mixture: {} with fractions {}".format(
            tissue_types, tissue_fractions
        )
    )

    fragments: FragmentList = simulator.simulate_tissue_mixture(
        tissue_types=tissue_types,
        tissue_fractions=tissue_fractions,
        total_fragments=total_fragments,
        genomic_regions=genomic_regions,
        add_noise=add_noise
    )

    logger.info(
        "saving {} simulated fragments to {}.frag".format(
            fragments.length(), output_name
        )
    )
    fragments.to_frag_file(output_name, out_dir)


def simulate_cancer_progression(
    config: dict[str, object], output_name: str, fasta_path: str, out_dir: str
) -> None:
    """Cancer progression simulation with multiple time points."""
    try:
        normal_profile: str = config["normal_profile"]  # type: ignore
        tumor_fractions: list[float] = \
            config["tumor_fractions"]  # type: ignore
        time_points: list[str] = config["time_points"]  # type: ignore
        fragments_per_timepoint: int = config[
            "fragments_per_timepoint"
        ]  # type: ignore
        bed_file: str = config["bed_file"]  # type: ignore
        if not os.path.exists(bed_file):
            fail("BED file '{}' does not exist".format(bed_file), logger)

        genomic_regions: list[tuple[str, int, int]] = parse_bed_file(bed_file)

    except KeyError as e:
        fail(
            "missing required parameter for cancer progression: {}".format(e),
            logger
        )

    logger.info(
        "initializing cancer progression simulator with FASTA: {}".format(
            fasta_path
        )
    )
    seq_gen = SequenceContextGenerator(fasta_path)
    simulator: TissueMixtureSimulator = TissueMixtureSimulator(seq_gen)

    logger.info(
        "simulating cancer progression: {} timepoints".format(
            len(time_points)
        )
    )

    results: dict[str, FragmentList] = (
        simulator.simulate_cancer_progression(
            normal_profile=normal_profile,
            tumor_fractions=tumor_fractions,
            time_points=time_points,
            fragments_per_timepoint=fragments_per_timepoint,
            genomic_regions=genomic_regions
        )
    )

    # @NOTE(ds): We save each time point as a separate .frag file.
    for timepoint, fragments in results.items():
        timepoint_name = "{}_{}".format(output_name, timepoint)
        logger.info(
            "saving {} fragments for timepoint {} to {}.frag".format(
                fragments.length(), timepoint, timepoint_name
            )
        )
        fragments.to_frag_file(timepoint_name, out_dir)


def switch_on_subcommand(subcmd: str, args: argparse.Namespace,
                         logger: logging.Logger) -> None:
    """
    Dispatch to appropriate subcommand handler.

    Routes execution to the correct function based on the subcommand.
    Handles output directory creation and logging setup.

    Args:
        subcmd: Subcommand name (extract, stats, lengths, scores, simulate,
            version)
        args: Parsed command line arguments
        logger: Logger instance for output

    Raises:
        CodeUnreachableError: If unknown subcommand is provided
    """
    if subcmd == "version":
        logger.info(version_string)
        exit(0)

    out_dir: str = args.out_dir
    if not os.path.isdir(out_dir):
        logger.info("creating directory `{}'".format(out_dir))
        os.makedirs(out_dir, exist_ok=True)
    logger.info("saving all output to `{}'".format(out_dir))

    if subcmd == "extract":
        extract(out_dir, args)
    elif subcmd == "stats":
        stats(out_dir, args)
    elif subcmd == "lengths":
        lengths(out_dir, args)
    elif subcmd == "scores":
        scores(out_dir, args)
    elif subcmd == "simulate":
        simulate(out_dir, args)
    else:
        raise CodeUnreachableError("unkown subcommand `{}'".format(subcmd))
    logger.info("done with `{}'".format(subcmd))


def search_dir(directory: str, exts: list[str]) -> list[str]:
    """
    Search for files with specific extensions in a directory.

    Args:
        directory: Directory path to search
        exts: List of file extensions to match (e.g., ['.bam', '.BAM'])

    Returns:
        List of full paths to matching files
    """
    results: list[str] = []
    potential_file: os.DirEntry[str]

    for potential_file in os.scandir(directory):
        if potential_file.is_file():
            ext: str
            file: str = potential_file.name
            _name, ext = os.path.splitext(file)
            if ext in exts:
                full_path: str = os.path.join(directory, file)
                results.append(full_path)

    return results


def get_frag_files(
    directory: Optional[str], file: Optional[str]
) -> list[str]:
    """
    Get list of fragment files from directory or single file.

    Args:
        directory: Directory containing .frag files (optional)
        file: Single .frag file path (optional)

    Returns:
        List of .frag file paths

    Raises:
        CodeUnreachableError: If neither directory nor file provided
        SystemExit: If directory/file doesn't exist
    """
    frag_files: list[str] = []
    if directory:
        if not os.path.isdir(directory):
            fail("directory `{}' does not exist".format(directory), logger)
        logger.info("looking for FRAG files in `{}'".format(directory))
        frag_files = search_dir(directory, [".FRAG", ".frag"])
    elif file:
        if not os.path.isfile(file):
            fail("file `{}' does not exist".format(file), logger)
        frag_files.append(file)
    else:
        raise CodeUnreachableError("internal error (no frag file or dir)")
    return frag_files


def filename_only(path: str) -> str:
    """
    Extract filename without extension from full path.

    Args:
        path: Full file path

    Returns:
        Filename without extension

    Note:
        Behavior is unspecified if path points to a directory
    """
    name: str
    basename: str = os.path.basename(path)
    name, _ = os.path.splitext(basename)
    return name


if __name__ == "__main__":
    logger: logging.Logger = logging.getLogger("pyfrag")
    argparser: argparse.ArgumentParser = create_argparser()
    args: argparse.Namespace = argparser.parse_args()

    level: int
    is_verbose: Final[bool] = args.is_verbose
    if is_verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.getLogger("pyfraglib").setLevel(level)
    logging.getLogger("pyfrag").setLevel(level)
    logging.getLogger("py.warnings").setLevel(level)

    signal.signal(signal.SIGINT, signal_handler)

    subcmd: str = args.subcommand
    switch_on_subcommand(subcmd, args, logger)
