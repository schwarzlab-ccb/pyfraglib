#!/usr/bin/env python3
"""
Transcription Start Site (TSS) Annotation Downloader

This script downloads transcription start site (TSS) annotations for a list of
genes from the Ensembl REST API and creates a BED file with TSS regions
±1000bp.

The script queries the Ensembl GRCh37/GRCh38 REST API to find canonical
transcripts for each gene and extracts the TSS coordinate based on the strand
orientation. For positive strand genes, TSS is the transcript start; for
negative strand genes, TSS is the transcript end.

Output format:
    BED file with columns: chromosome, start, end, gene_symbol
    Regions are TSS ±1000bp, sorted by chromosome order

Usage:
    python download_tss_annos.py -g genes.txt -o output.bed --ref hg19

Input file format:
    Text file with one gene symbol per line (e.g., EGFR, TP53, BRCA1)

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
import requests
import logging
import os
import re
import sys
import tqdm
import pyfraglib

from typing import Final, NoReturn

version_string: Final[str] = "download_tss_annos v{} (running on Python v{})" \
    .format(pyfraglib.__version__, sys.version.split(" ")[0])
LOGGER_NAME: Final[str] = "download_tss_annos"


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
    Create and configure the argument parser for the TSS downloader.

    Returns:
        Configured ArgumentParser instance with required options:
        - gene_list: Input file with gene symbols
        - ref: Genome assembly (hg19 or hg38, default: hg19)
        - outfile: Output BED file path
        - verbose: Enable debug logging
    """
    argparser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="download_tss_annos", description="Download transcription "
        "start site (TSS) annotations for a list of genes.",
        epilog="{}. Licensed under GPLv3. See repository at `{}' for more "
        "info.".format(version_string, pyfraglib.__repo_url__))
    argparser.add_argument(
        "-g", "--gene-list", type=str, dest="genes_file", required=True,
        help="Gene list TXT file (one gene name per line).")
    argparser.add_argument(
        "--ref", type=str, dest="genome_assembly", required=False,
        default="hg19", help="Genome assembly for which to download TSSs.")
    argparser.add_argument(
        "-o", "--outfile", type=str, dest="outfile", required=True,
        help="The name of the BED output file.")
    argparser.add_argument(
        "-v", "--verbose", action="store_true", dest="is_verbose",
        default=False, help="Set log level to print debugging info.")

    return argparser


def get_canonical_tss(
    gene_symbol: str, assembly: str = "hg19"
) -> dict[str, object]:
    """
    Retrieve canonical transcript TSS information for a gene from Ensembl.

    Queries the Ensembl GRCh37/GRCh38 REST API to find the canonical transcript
    for the given gene symbol and extracts TSS coordinates.

    Args:
        gene_symbol: Gene symbol (e.g., 'EGFR', 'TP53')
        assembly: Genome assembly (either hg19 or hg38)

    Returns:
        Dictionary containing:
        - gene: Gene symbol
        - start: Transcript start coordinate
        - end: Transcript end coordinate
        - canonical_transcript: Ensembl transcript ID
        - contig: Chromosome name
        - strand: Strand orientation (1 or -1)
        - TSS: Transcription start site coordinate

    Raises:
        ValueError: If no canonical transcript found for the gene

    Note:
        TSS is transcript start for positive strand, transcript end for
        negative strand
    """
    api_base_hg19: Final[str] = "https://grch37.rest.ensembl.org"
    api_base_hg38: Final[str] = "https://rest.ensembl.org"

    api_base: str
    if assembly == "hg19":
        api_base = api_base_hg19
    elif assembly == "hg38":
        api_base = api_base_hg38
    else:
        raise ValueError(f"unknown genome assembly {assembly}")

    url: str = f"{api_base}/lookup/symbol/homo_sapiens/" \
               f"{gene_symbol}?expand=1"
    headers: dict[str, str] = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)
    data: object = response.json()

    transcript: dict[str, object]
    for transcript in data.get("Transcript", []):  # type: ignore
        if transcript.get("is_canonical"):
            return {
                "gene": gene_symbol,
                "start": transcript["start"],
                "end": transcript["end"],
                "canonical_transcript": transcript["id"],
                "contig": transcript["seq_region_name"],
                "strand": transcript["strand"],
                "TSS": transcript["start"]
                if transcript["strand"] == 1 else transcript["end"]
            }
    raise ValueError(f"no canonical transcript found for {gene_symbol}")


def chrom_order(name: str) -> int:
    """
    Get chromosome sorting order for BED file output.

    Converts chromosome names to numeric values for proper sorting:
    - Chromosomes 1-22: return numeric value
    - X chromosome: return 23
    - Y chromosome: return 24
    - M chromosome: return 25
    - Other chromosomes: return 0

    Args:
        name: Chromosome name (e.g., 'chr1', '1', 'chrX', 'X')

    Returns:
        Numeric value for chromosome sorting

    Note:
        Handles both 'chr' prefixed and non-prefixed chromosome names
    """
    match = re.match(r"(?:chr)?(\d+|X|Y)", name)
    if not match:
        return 0

    val: str = match.group(1)
    if val == "X":
        return 23
    elif val == "Y":
        return 24
    elif val == "M":
        return 25
    else:
        return int(val)


if __name__ == "__main__":
    """
    Main execution block for TSS annotation downloader.

    Processes command line arguments, reads gene list, queries Ensembl API
    for each gene, and creates a sorted BED file with TSS regions.

    Workflow:
    1. Parse command line arguments
    2. Read gene list from input file
    3. Query Ensembl API for each gene's canonical TSS
    4. Create BED records with TSS ±1000bp regions
    5. Sort records by chromosome order
    6. Write sorted BED file

    Exit codes:
    - 0: Success
    - 1: Error (file not found, gene lookup failed, etc.)
    """
    logger: logging.Logger = logging.getLogger(LOGGER_NAME)
    argparser: argparse.ArgumentParser = create_argparser()
    args: argparse.Namespace = argparser.parse_args()

    level: int
    is_verbose: Final[bool] = args.is_verbose
    if is_verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logger.setLevel(level)
    logging.getLogger("urllib3.connectionpool").setLevel(level)

    genes_file: Final[str] = args.genes_file
    outfile: Final[str] = args.outfile
    genome_assembly: Final[str] = args.genome_assembly

    if not os.path.isfile(genes_file):
        fail("gene list `{}' does not exist".format(genes_file), logger)
    elif os.path.isfile(outfile):
        fail("file `{}' already exists".format(outfile), logger)
    elif genome_assembly not in ["hg19", "hg38"]:
        fail("genome assembly must be 'hg19' or 'hg38', got '{}'".format(
            genome_assembly), logger)

    gene_list: list[str]
    with open(genes_file, "r") as genes:
        gene_list = [gene.strip("\n") for gene in genes.readlines()]

    bed_records: list[tuple[object, int, int, object]] = []
    for gene in tqdm.tqdm(gene_list, desc="Downloading TSS annotations"):
        try:
            d: dict[str, object] = get_canonical_tss(gene, genome_assembly)
            tss: int = int(d["TSS"])  # type: ignore
            bed_records.append(
                (d["contig"], tss-1000, tss+1000, d["gene"])
            )
            logger.debug(f"TSS found for {gene}")
        except ValueError as e:
            logger.warning(e)

    # @NOTE(ds): Sort BED records by chromosome order, then by position.
    bed_records.sort(
        key=lambda x: (chrom_order(x[0]), x[1])  # type: ignore
    )
    logger.info(f"Writing {len(bed_records)} TSS regions to {outfile}")
    with open(outfile, "w") as bed_file:
        for record in bed_records:
            bed_file.write(
                f"{record[0]}\t{record[1]}\t{record[2]}\t{record[3]}\n"
            )
    logger.info("TSS annotation download completed successfully")
