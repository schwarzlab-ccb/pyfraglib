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
import codecs
import os

from setuptools import setup
from typing import cast


def read(rel_path: str) -> str:
    here: str = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path: str) -> str:
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name="pyfraglib",
    version=get_version("pyfraglib/__init__.py"),
    url="https://bitbucket.org/schwarzlab/project-lymphoma-cfdna/",
    author="Daniel Schütte",
    author_email="daniel.schuette@iccb-cologne.org",
    description="Software suite to calculate fragmentomics features "
                "from cfDNA and perform downstream analyses.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    license="GPL-v3",
    packages=["pyfraglib", "pyfraglib.simulator"],
    install_requires=cast(list[str], [
        "pysam==0.22.1", "matplotlib==3.9.0", "numpy<1.27", "seaborn==0.13.2",
        "matplotlib-stubs==0.2.0", "types-seaborn==0.13.2.20240820",
        "pandas==2.2.2", "intervaltree==3.1.0", "scipy==1.14.1",
        "scipy-stubs==1.14.1.4", "tqdm==4.67.0", "types-tqdm==4.67.0.20241119",
        "mypy==1.13.0", "flake8==7.1.1"
    ]),
    package_data={"pyfraglib": ["py.typed"]},
    scripts=[
        "scripts/pyfrag.py",
        "scripts/txt_to_vcf.py",
        "scripts/download_tss_annos.py",
        "scripts/extract_mutated_reads.py"
    ],

    keywords=["cfDNA", "fragmentomics"],
    classifiers=[
        "Programming Language :: Python :: 3"
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
