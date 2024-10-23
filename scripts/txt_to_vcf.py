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
import re
import sys
import pyfraglib
import pysam

from typing import Final, NoReturn
from pyfraglib.core import hg19_chromosomes

version_string: Final[str] = "txt_to_vcf v{} (running on Python v{})" \
    .format(pyfraglib.__version__, sys.version.split(" ")[0])


# @NOTE(ds): We re-define fail to be more specific with the logger that we
# use. Otherwise, we could have just used pyfraglib's `fail'.
def fail(msg: str, logger: logging.Logger) -> NoReturn:
    logger.fatal(msg)
    sys.exit(1)


def create_argparser() -> argparse.ArgumentParser:
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
        "-v", "--verbose", action="store_true", dest="is_verbose",
        default=False, help="Set log level so that debugging info is printed.")

    return argparser


def convert(infile: str, outfile: str, logger: logging.Logger) -> None:
    vcf_header: pysam.VariantHeader = pysam.VariantHeader()
    vcf_header.add_meta("fileformat", "VCFv4.2")
    vcf_header.add_meta("source", "pyfraglib")

    for chrom, length, _, _ in hg19_chromosomes:
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
                    # invalid.
                    info_name: str = field_names[it].replace("-", "_")
                    info_name = re.sub("[+-]", "*", info_name)
                    field_names[it] = info_name

                    vcf_file.header.info.add(
                        info_name, number=".", type="String", description=""
                    )
            else:
                new_record: pysam.VariantRecord = vcf_file.new_record()

                new_record.chrom = fields[0]
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

    if not os.path.isfile(infile):
        fail("file `{}' does not exist".format(infile), logger)
    elif os.path.isfile(outfile):
        fail("file `{}' already exists".format(outfile), logger)

    convert(infile, outfile, logger)
