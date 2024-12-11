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
import logging
import os
import signal
import sys

import numpy as np

from multiprocessing.managers import BaseManager
from typing import NoReturn, Final

# @NOTE(ds): From ``https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37''.
hg19_chromosomes: Final[list[tuple[str, int, str, str]]] = [
    ("1",  249250621, "CM000663.1", "NC_000001.10"),
    ("2",  243199373, "CM000664.1", "NC_000002.11"),
    ("3",  198022430, "CM000665.1", "NC_000003.11"),
    ("4",  191154276, "CM000666.1", "NC_000004.11"),
    ("5",  180915260, "CM000667.1", "NC_000005.9"),
    ("6",  171115067, "CM000668.1", "NC_000006.11"),
    ("7",  159138663, "CM000669.1", "NC_000007.13"),
    ("8",  146364022, "CM000670.1", "NC_000008.10"),
    ("9",  141213431, "CM000671.1", "NC_000009.11"),
    ("10", 135534747, "CM000672.1", "NC_000010.10"),
    ("11", 135006516, "CM000673.1", "NC_000011.9"),
    ("12", 133851895, "CM000674.1", "NC_000012.11"),
    ("13", 115169878, "CM000675.1", "NC_000013.10"),
    ("14", 107349540, "CM000676.1", "NC_000014.8"),
    ("15", 102531392, "CM000677.1", "NC_000015.9"),
    ("16", 90354753,  "CM000678.1", "NC_000016.9"),
    ("17", 81195210,  "CM000679.1", "NC_000017.10"),
    ("18", 78077248,  "CM000680.1", "NC_000018.9"),
    ("19", 59128983,  "CM000681.1", "NC_000019.9"),
    ("20", 63025520,  "CM000682.1", "NC_000020.10"),
    ("21", 48129895,  "CM000683.1", "NC_000021.8"),
    ("22", 51304566,  "CM000684.1", "NC_000022.10"),
    ("X",  155270560, "CM000685.1", "NC_000023.10"),
    ("Y",  59373566,  "CM000686.1", "NC_000024.9"),
    ("M",  16569,     "J01415.2",   "NC_012920.1")
]

# @NOTE(ds): From ``https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38''.
hg38_chromosomes: Final[list[tuple[str, int, str, str]]] = [
    ("1",  248956422, "CM000663.2", "NC_000001.11"),
    ("2",  242193529, "CM000664.2", "NC_000002.12"),
    ("3",  198295559, "CM000665.2", "NC_000003.12"),
    ("4",  190214555, "CM000666.2", "NC_000004.12"),
    ("5",  181538259, "CM000667.2", "NC_000005.10"),
    ("6",  170805979, "CM000668.2", "NC_000006.12"),
    ("7",  159345973, "CM000669.2", "NC_000007.14"),
    ("8",  145138636, "CM000670.2", "NC_000008.11"),
    ("9",  138394717, "CM000671.2", "NC_000009.12"),
    ("10", 133797422, "CM000672.2", "NC_000010.11"),
    ("11", 135086622, "CM000673.2", "NC_000011.10"),
    ("12", 133275309, "CM000674.2", "NC_000012.12"),
    ("13", 114364328, "CM000675.2", "NC_000013.11"),
    ("14", 107043718, "CM000676.2", "NC_000014.9"),
    ("15", 101991189, "CM000677.2", "NC_000015.10"),
    ("16", 90338345,  "CM000678.2", "NC_000016.10"),
    ("17", 83257441,  "CM000679.2", "NC_000017.11"),
    ("18", 80373285,  "CM000680.2", "NC_000018.10"),
    ("19", 58617616,  "CM000681.2", "NC_000019.10"),
    ("20", 64444167,  "CM000682.2", "NC_000020.11"),
    ("21", 46709983,  "CM000683.2", "NC_000021.9"),
    ("22", 50818468,  "CM000684.2", "NC_000022.11"),
    ("X",  156040895, "CM000685.2", "NC_000023.11"),
    ("Y",  57227415,  "CM000686.2", "NC_000024.10"),
    ("M",  16569,     "J01415.2",   "NC_012920.1")
]


# `fail' terminates the program even if called from a subprocess.
def fail(msg: str) -> NoReturn:
    logger: logging.Logger = logging.getLogger("pyfraglib")
    logger.fatal(msg)

    # @NOTE(ds): We signal to a potential parent.
    os.kill(os.getppid(), signal.SIGINT)

    sys.exit(1)


class PyfraglibException(Exception):
    def __init__(self, msg: str) -> None:
        logger: logging.Logger = logging.getLogger("pyfraglib")
        logger.fatal(msg)


class CodeUnreachableError(PyfraglibException):
    def __init__(self, msg: str) -> None:
        super().__init__(msg)


class PyfragManager(BaseManager):
    pass


def get_chromosome_length(chrom: str, genome: str = "hg19") -> int:
    chromosomes: list[tuple[str, int, str, str]]
    if genome == "hg19":
        chromosomes = hg19_chromosomes
    elif genome == "hg38":
        chromosomes = hg38_chromosomes
    else:
        fail("unknown genome `{}' requested".format(genome))

    chrom = homogenize_contig_name(chrom)
    for name, length, _, _ in chromosomes:
        if chrom == name:
            return length
    fail("unknown chromosome `{}'".format(chrom))


def homogenize_contig_name(contig: str) -> str:
    if contig.startswith("chr"):
        contig = contig[3:]
    return contig


# @NOTE(ds): The input list needs to be normalized to proportions!
def shannon_entropy(proportions: list[float]) -> float:
    prop_sum: float = 0.0
    for prop in proportions:
        log_prop: float = np.log(prop)
        prop_sum += log_prop * prop

    # @NOTE(ds): Scipy can calculate entropy, too:
    #  > from scipy.stats import entropy
    #  >
    #  > scipy_shannon: float = typing.cast(float, entropy(pk=proportions))
    return -prop_sum


def simpson_index(proportions: list[float]) -> float:
    prop_sum: float = 0.0
    for prop in proportions:
        prop_sum += prop * prop
    return prop_sum


def detect_cpus() -> int:
    logger: logging.Logger = logging.getLogger("pyfraglib")
    num_cores: Final[int] = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))

    logger.info("{} core(s) detected".format(num_cores))

    return num_cores
