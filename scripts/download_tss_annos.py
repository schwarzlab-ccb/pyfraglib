#!/usr/bin/env python3
#
# This file is part of `pyfraglib`, a software suite to calculate fragmentomics
# features from cfDNA and perform downstream analyses.
#
# Copyright (C) 2025 Daniel Schütte, daniel.schuette@iccb-cologne.org
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


# @NOTE(ds): We re-define fail to be more specific with the logger that we
# use. Otherwise, we could have just used pyfraglib's `fail'.
def fail(msg: str, logger: logging.Logger) -> NoReturn:
    logger.fatal(msg)
    sys.exit(1)


def create_argparser() -> argparse.ArgumentParser:
    argparser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="download_tss_annos", description="Download transcription "
        "start site (TSS) annotations for a list of genes.",
        epilog="{}. Licensed under GPLv3. See repository at `{}' for more "
        "info.".format(version_string, pyfraglib.__repo_url__))
    argparser.add_argument(
        "-g", "--gene-list", type=str, dest="genes_file", required=True,
        help="Gene list TXT file (one gene name per line).")
    argparser.add_argument(
        "-o", "--outfile", type=str, dest="outfile", required=True,
        help="The name of the BED output file.")
    argparser.add_argument(
        "-v", "--verbose", action="store_true", dest="is_verbose",
        default=False, help="Set log level so that debugging info is printed.")

    return argparser


def get_canonical_tss(gene_symbol: str) -> dict[str, object]:
    url: str = "https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/" \
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
    logger: logging.Logger = logging.getLogger("download_tss_annos")
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

    if not os.path.isfile(genes_file):
        fail("gene list `{}' does not exist".format(genes_file), logger)
    elif os.path.isfile(outfile):
        fail("file `{}' already exists".format(outfile), logger)

    gene_list: list[str]
    with open(genes_file, "r") as genes:
        gene_list = [gene.strip("\n") for gene in genes.readlines()]

    bed_records: list[tuple[object, int, int, object]] = []
    for gene in tqdm.tqdm(gene_list):
        try:
            d: dict[str, object] = get_canonical_tss(gene)
            tss: int = int(d["TSS"])  # type: ignore
            bed_records.append(
                (d["contig"], tss-1000, tss+1000, d["gene"])
            )
            logger.debug(f"TSS found for {gene}")
        except ValueError as e:
            logger.warning(e)

    bed_records.sort(
        key=lambda x: (chrom_order(x[0]), x[1])  # type: ignore
    )
    with open(outfile, "w") as bed_file:
        for record in bed_records:
            bed_file.write(
                f"{record[0]}\t{record[1]}\t{record[2]}\t{record[3]}\n"
            )
