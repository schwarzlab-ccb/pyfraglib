#!/usr/bin/env python3
"""
Extract mutated reads from a BAM file based on VCF variant information.

This script subsets a BAM file to contain only reads that support variants
defined in a VCF file. It uses the same logic as the Fragment class for
consistent mutation annotation.

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
import pysam
import pyfraglib

from tqdm import tqdm
from typing import Final, NoReturn, Set
from pyfraglib.core import get_logger, homogenize_to_chrom_naming_convention

version_string: Final[str] = \
    "extract_mutated_reads v{} (running on Python v{})" \
    .format(pyfraglib.__version__, sys.version.split(" ")[0])
LOGGER_NAME: Final[str] = "extract_mutated_reads"


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


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        prog="extract_mutated_reads", description="Extract mutated reads from "
        "BAM file based on VCF variants.",
        epilog="{}. Licensed under GPLv3. See repository at `{}' for more "
        "info.".format(version_string, pyfraglib.__repo_url__)
    )
    parser.add_argument(
        "-b", "--bam", required=True, help="Input BAM file path."
    )
    parser.add_argument(
        "-m", "--vcf", required=True, help="Input VCF file path."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output BAM file path."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging."
    )

    return parser.parse_args()


def load_variants(
    vcf_path: str, logger: logging.Logger
) -> dict[tuple[str, int], set[str]]:
    """
    Load variants from VCF file.

    Args:
        vcf_path: Path to VCF file
        logger: Logger to use

    Returns:
        Dictionary mapping (chromosome, position) to set of alternative alleles
    """
    logger = get_logger()
    variants: dict[tuple[str, int], set[str]] = {}

    try:
        with pysam.VariantFile(vcf_path) as vcf:
            rec: pysam.VariantRecord
            for rec in vcf.fetch():
                if rec.rlen != 1:  # only process SNVs
                    continue
                if rec.alts is not None:
                    # @NOTE(ds): `pysam' has 2 APIs for accessing the genomic
                    # position of a variant: whereas `pos' is 1-based, `start'
                    # is 0-based. Thus, we must be very careful here. See also:
                    # ``https://pysam.readthedocs.io/en/latest/''.
                    variants.setdefault(
                        (rec.contig, rec.start), set()
                    ).update(rec.alts)
    except Exception as e:
        fail(f"Error reading VCF file '{vcf_path}': {e}", logger)

    logger.info(f"Loaded {len(variants)} variants from {vcf_path}")
    return variants


def find_supporting_reads(
    bam_path: str, variants: dict[tuple[str, int], set[str]],
    bam_header: dict[str, object], logger: logging.Logger
) -> set[str]:
    """
    Find reads that support variants.

    Args:
        bam_path: Path to BAM file
        variants: Dictionary of variants from VCF
        bam_header: BAM file header for chromosome normalization
        logger: Logger to use

    Returns:
        Set of read names that support variants
    """
    supporting_reads: set[str] = set()
    normalized_variants: dict[tuple[str, int], set[str]] = {}
    for (chrom, pos), alts in variants.items():
        normalized_chrom = homogenize_to_chrom_naming_convention(
            chrom, bam_header
        )
        normalized_variants[(normalized_chrom, pos)] = alts

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        total_reads: int = 0
        processed_reads: int = 0

        for read in tqdm(bam.fetch(until_eof=True), desc="Processing reads"):
            total_reads += 1

            if (read.is_unmapped or read.query_name is None or
                    read.query_sequence is None):
                continue

            processed_reads += 1
            chrom = bam.get_reference_name(read.reference_id)
            if chrom is None:
                continue  # type: ignore[unreachable]

            read_positions: list[int | None] = \
                read.get_reference_positions(full_length=True)  # type: ignore

            for qpos, ref_pos in enumerate(read_positions):
                if ref_pos is None:
                    continue  # skip this position bc. it's a gap/deletion

                key = (chrom, ref_pos)  # 0-based position
                if key in normalized_variants:
                    read_base = read.query_sequence[qpos]
                    if read_base in normalized_variants[key]:
                        supporting_reads.add(read.query_name)
                        break

    logger.info(f"Processed {processed_reads} valid reads out of "
                f"{total_reads} total reads")
    logger.info(f"Found {len(supporting_reads)} reads supporting variants")
    return supporting_reads


def extract_mutated_reads(
    bam_path: str, output_path: str, supporting_reads: Set[str],
    logger: logging.Logger
) -> None:
    """
    Extract supporting reads to output BAM file.

    Args:
        bam_path: Input BAM file path
        output_path: Output BAM file path
        supporting_reads: Set of read names to extract
        logger: Logger to use
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        with pysam.AlignmentFile(output_path, "wb", template=bam) as output:
            total_reads: int = 0
            extracted_reads: int = 0

            for read in tqdm(
                bam.fetch(until_eof=True), desc="Extracting reads"
            ):
                total_reads += 1
                if read.query_name in supporting_reads:
                    output.write(read)
                    extracted_reads += 1

    logger.info(f"Extracted {extracted_reads} reads out of "
                f"{total_reads} total reads")
    logger.info(f"Output written to {output_path}")


def main() -> None:
    """Main function."""
    logger: logging.Logger = logging.getLogger(LOGGER_NAME)
    args: argparse.Namespace = parse_arguments()

    bam_path: str = args.bam
    vcf_path: str = args.vcf
    output_path: str = args.output
    verbose: bool = args.verbose

    if verbose:
        logger.setLevel(logging.DEBUG)

    if not os.path.exists(bam_path):
        fail(f"BAM file not found: {bam_path}", logger)
    elif not os.path.exists(vcf_path):
        fail(f"VCF file not found: {vcf_path}", logger)
    elif os.path.exists(output_path):
        fail(f"Output file already exists: {output_path}", logger)

    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            if not bam.has_index():
                fail(f"BAM file '{bam_path}' is not indexed. "
                     f"Please create an index file.", logger)
            bam_header: dict[str, object] = bam.header.to_dict()
    except Exception as e:
        fail(f"Error reading BAM file '{bam_path}': {e}", logger)

    output_dir: str = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    variants = load_variants(vcf_path, logger)
    if not variants:
        fail("No variants found in VCF file", logger)

    supporting_reads = find_supporting_reads(
        bam_path, variants, bam_header, logger
    )
    if not supporting_reads:
        fail("No reads found supporting variants", logger)

    extract_mutated_reads(bam_path, output_path, supporting_reads, logger)
    logger.info("Processing complete")


if __name__ == "__main__":
    main()
