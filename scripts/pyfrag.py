#!/usr/bin/env python3
#
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
import argparse
import logging
import os
import signal
import sys
import pysam

import pandas as pd

from typing import Final, NoReturn, Optional

import pyfraglib
from pyfraglib import Fragment, FragmentList
from pyfraglib.core import CodeUnreachableError
from pyfraglib.fragfile import FragFile
from pyfraglib.lengths import fragment_length_plot, fragment_length_gmm
from pyfraglib.stats import fragments_per_chromosome_barplot, \
                            end_motifs_barplot, log_stats
from pyfraglib.scores import motif_diversity, windowed_protection_score, \
                             score_line_plot

version_string: Final[str] = "pyfraglib v{} (running on Python v{})" \
    .format(pyfraglib.__version__, sys.version.split(" ")[0])


# @NOTE(ds): We re-define fail to be more specific with the logger that we
# use. Otherwise, we could have just used pyfraglib's `fail'. Also, we
# probably do not need to signal a parent process.
def fail(msg: str, logger: logging.Logger) -> NoReturn:
    logger.fatal(msg)
    sys.exit(1)


def signal_handler(sig: int, frame: object) -> NoReturn:
    logger: logging.Logger = logging.getLogger("pyfrag")
    fail("an error occurred", logger)


# @NOTE(ds): We always require our users to indicate an output directory. What
# type of output we produce is determined by a subcommand (implemented via
# argparse's subparsers). Subcommands have a bunch of additional, but
# individual options which determine their behavior only.
def create_argparser() -> argparse.ArgumentParser:
    argparser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="pyfrag", description="Use a subset of `pyfraglib's "
        "capabilities from the command line.",
        epilog="{}. Licensed under GPLv3. See repository at `{}' for more "
        "info.".format(version_string, pyfraglib.__repo_url__))
    argparser.add_argument(
        "-o", "--out-dir", type=str, dest="out_dir", default="pyfrag_out",
        help="Directory that output is saved to. The type of output depends "
        "on the subcommand. Dir is created if it does not exist.")
    argparser.add_argument(
        "-v", "--verbose", action="store_true", dest="is_verbose",
        default=False, help="Set log level so that debugging info is printed.")

    subparsers: argparse._SubParsersAction[argparse.ArgumentParser] = \
        argparser.add_subparsers(
            required=True, help="Select from a number of utilities. "
            "Example: ``$ pyfraglib.py extract --bam-file=<FILE>''",
            dest="subcommand")

    subparsers.add_parser("version", help="Show version info.")

    argparser_extract: argparse.ArgumentParser = subparsers.add_parser(
        "extract", help="Extract fragment information from BAM file(s). "
        "Results are written to disk in the `.frag' file format.")
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
        "statistics and save them to disk.")
    argparser_stats.add_argument(
        "-d", "--frag-dir", type=str, dest="frag_dir", required=False,
        help="Input FRAG files to be analyzed in batch mode. Must not be "
        "combined with `--frag-file'. Requires sufficient amounts of "
        "memory.")
    argparser_stats.add_argument(
        "-f", "--frag-file", type=str, dest="frag_file", required=False,
        help="Input FRAG file to be analyzed. Must not be combined with "
        "`--frag-dir'.")

    argparser_lengths: argparse.ArgumentParser = subparsers.add_parser(
        "lengths", help="Given (a) `.frag' file(s), read fragment length "
        "data and create a number of plots that are saved to disk.")
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
        "fragmentomics scores and save the results to disk.")
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

    return argparser


# @NOTE(ds): At this point, we can assume that `out_dir' exists.
# Loading all `FragmentList's into a collection before saving the to disk is
# super expensive in terms of memory. We did that in previous versions of
# `pyfraglib' but changed things around to be able to handle large directories
# of even larger BAM files.
def extract(out_dir: str, args: argparse.Namespace) -> None:
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


# @NOTE(ds): `stats' and `lengths' follow very similar patterns.
def stats(out_dir: str, args: argparse.Namespace) -> None:
    frag_file: Final[str] = args.frag_file
    frag_dir: Final[str] = args.frag_dir

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
        end_motifs_barplot(fragments, out_dir, name, kmer_len=3)


def lengths(out_dir: str, args: argparse.Namespace) -> None:
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


# @NOTE(ds): We interleave the calculation of our scores. It's ugly, but
# otherwise we have to store all fragment files as a collection which requires
# loads of memory.
def scores(out_dir: str, args: argparse.Namespace) -> None:
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


def switch_on_subcommand(subcmd: str, args: argparse.Namespace,
                         logger: logging.Logger) -> None:
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
    else:
        raise CodeUnreachableError("unkown subcommand `{}'".format(subcmd))
    logger.info("done with `{}'".format(subcmd))


def search_dir(directory: str, exts: list[str]) -> list[str]:
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


# @NOTE(ds): Given a full path, extract only the name of file without its
# extension. If the last component of the path is _not_ a file but a directory,
# the behavior of this function is unspecified.
def filename_only(path: str) -> str:
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
