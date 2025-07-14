#!/usr/bin/env python3
"""
TXT to VCF Converter

This script converts custom tab-delimited SNV files to standard VCF format.
It processes files with genomic variant information and creates properly
formatted VCF files compatible with standard bioinformatics tools.

Input format:
    Tab-delimited file with required columns:
    - Chr: Chromosome name
    - Start: Start position (1-based)
    - End: End position (1-based)
    - Ref: Reference allele
    - Alt: Alternative allele

    Additional columns are preserved as INFO fields in the VCF output.
    The second line may contain column descriptions that will be used
    as INFO field names.

Output format:
    Standard VCF 4.2 format with:
    - Proper VCF header with contigs for specified reference genome
    - All additional fields stored as INFO annotations
    - FILTER field set to PASS for all variants

Usage:
    python txt_to_vcf.py -f input.txt -o output.vcf -g hg19

Supported reference genomes:
    - hg19: Human genome build 19 (GRCh37)
    - hg38: Human genome build 38 (GRCh38)

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
import re
import sys
import pyfraglib
import pysam

from typing import Final, NoReturn
from pyfraglib.core import hg19_chromosomes, hg38_chromosomes, \
                           homogenize_contig_name

version_string: Final[str] = "txt_to_vcf v{} (running on Python v{})" \
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


def create_argparser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser for the TXT to VCF converter.

    Returns:
        Configured ArgumentParser instance with required options:
        - infile: Input TXT file path
        - outfile: Output VCF file path
        - ref_genome: Reference genome (hg19 or hg38)
        - verbose: Enable debug logging
    """
    argparser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="txt_to_vcf", description="Convert our custom SNV file format to "
        "VCF format.",
        epilog="{}. Licensed under GPLv3. See repository at `{}' for more "
        "info.".format(version_string, pyfraglib.__repo_url__))
    argparser.add_argument(
        "-f", "--infile", type=str, dest="infile", required=True,
        help="The name of the TXT input file.")
    argparser.add_argument(
        "-o", "--outfile", type=str, dest="outfile", required=True,
        help="The name of the VCF output file.")
    argparser.add_argument(
        "-g", "--ref-genome", type=str, dest="ref_genome", required=True,
        help="Indicate the ref genome you are using (must be hg19 or hg38).")
    argparser.add_argument(
        "-v", "--verbose", action="store_true", dest="is_verbose",
        default=False, help="Set log level so that debugging info is printed.")

    return argparser


def convert(
    infile: str, outfile: str, ref_genome: str, logger: logging.Logger
) -> None:
    """
    Convert custom TXT format to standard VCF format.

    Processes a tab-delimited file with genomic variants and creates a
    properly formatted VCF file. Handles column mapping, INFO field
    creation, and contig definitions based on reference genome.

    Args:
        infile: Path to input TXT file
        outfile: Path to output VCF file
        ref_genome: Reference genome version ('hg19' or 'hg38')
        logger: Logger instance for status messages

    File format requirements:
        - Line 1: Header with column names (Chr, Start, End, Ref, Alt required)
        - Line 2: Optional column descriptions (used as INFO field names)
        - Line 3+: Variant data rows

    Raises:
        SystemExit: If required columns are missing or file format is invalid

    Note:
        - Chromosome names are normalized using homogenize_contig_name
        - All variants are marked as PASS in the FILTER field
        - Additional columns become INFO fields with String type
        - Invalid characters in field names are replaced with placeholders
    """
    vcf_header: pysam.VariantHeader = pysam.VariantHeader()
    vcf_header.add_meta("fileformat", "VCFv4.2")
    vcf_header.add_meta("source", "pyfraglib")

    genome: list[tuple[str, int, str, str]] = \
        hg19_chromosomes if ref_genome == "hg19" else hg38_chromosomes
    for chrom, length, _, _ in genome:
        vcf_header.contigs.add(chrom, length=length)

    vcf_file: pysam.VariantFile = \
        pysam.VariantFile(outfile, "w", header=vcf_header)

    field_names: list[str]
    required_fields: list[str] = ["Chr", "Start", "End", "Ref", "Alt"]
    with open(infile, "r") as txt_file:
        line: str
        for it, line in enumerate(txt_file.readlines()):
            line = line.rstrip("\n")
            fields: list[str] = line.split("\t")

            if it == 0:
                field_names = fields
                if not all([fn in field_names for fn in required_fields]):
                    fail("missing some required fields ({})".format(
                        required_fields), logger)

                # Also, we make some assumptions with regards to field naming
                # and column contents.
                assert field_names.index("Chr") == 0
                assert field_names.index("Start") == 1
                assert field_names.index("End") == 2
                assert field_names.index("Ref") == 3
                assert field_names.index("Alt") == 4

            elif it == 1:
                # This is a weird second header. We substitute its meaningful
                # labels into our field names.
                for it, name in enumerate(fields):
                    if it < len(required_fields):
                        continue

                    if name.upper() == "NA" or name == "":
                        pass
                    else:
                        field_names[it] = name

                    # We know very little about the info that our fields
                    # provide. And unfortunately, some of the info names are
                    # invalid. We substitute the invalid parts with a
                    # placeholder to avoid warnings later on.
                    info_name: str = field_names[it].replace("-", "_")
                    info_name = re.sub("[*+-]", "INVALID{}".format(it),
                                       info_name)
                    field_names[it] = info_name

                    vcf_file.header.info.add(
                        info_name, number=".", type="String", description=""
                    )
            else:
                new_record: pysam.VariantRecord = vcf_file.new_record()

                new_record.chrom = homogenize_contig_name(fields[0])
                new_record.pos = int(fields[1])  # indexing!
                new_record.rlen = int(fields[2]) - new_record.pos + 1
                new_record.ref = fields[3]
                new_record.filter.add("PASS")

                assert len(fields[4]) == 1  # always just 1 alt allele
                new_record.alts = (fields[4],)

                # Add all other fields to INFO. Some fields like VAF could
                # probably go into a more specific VCF column, but for now we
                # don't care.
                for it, field_name in enumerate(field_names):
                    if it < len(required_fields):
                        continue
                    info: str = fields[it].replace(";", ":")
                    new_record.info[field_name] = info

                vcf_file.write(new_record)


if __name__ == "__main__":
    """
    Main execution block for TXT to VCF converter.

    Processes command line arguments, validates input/output files,
    and performs the conversion from custom TXT format to standard VCF.

    Workflow:
    1. Parse command line arguments
    2. Validate input file exists and output file doesn't exist
    3. Call convert function to perform the conversion
    4. Log completion status

    Exit codes:
    - 0: Success
    - 1: Error (file not found, invalid format, etc.)
    """
    logger: logging.Logger = logging.getLogger("txt_to_vcf")
    argparser: argparse.ArgumentParser = create_argparser()
    args: argparse.Namespace = argparser.parse_args()

    level: int
    is_verbose: Final[bool] = args.is_verbose
    if is_verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logger.setLevel(level)

    infile: Final[str] = args.infile
    outfile: Final[str] = args.outfile
    ref_genome: Final[str] = args.ref_genome

    if not os.path.isfile(infile):
        fail("file `{}' does not exist".format(infile), logger)
    elif os.path.isfile(outfile):
        fail("file `{}' already exists".format(outfile), logger)
    elif ref_genome not in ["hg19", "hg38"]:
        fail("reference genome must be 'hg19' or 'hg38', got '{}'".format(
            ref_genome), logger)

    logger.info("Converting {} to VCF format using {} reference".format(
        infile, ref_genome))
    convert(infile, outfile, ref_genome, logger)
    logger.info("Conversion completed successfully: {}".format(outfile))
