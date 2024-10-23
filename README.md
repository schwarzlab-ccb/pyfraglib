# PyFragLib

## Overview

`pyfraglib` is a Python library to analyze high-throughput sequencing data of
cell-free DNA (cfDNA). More specifically, it facilitates the investigation of
fragmentomics features of which a list can be found below.

Because fragmentation of cfDNA is non-random, tissue- and disease-specific
patterns emerge in e.g. fragment length or end motif distributions. To date,
no comprehensive tool exists to implement a "fragmentomics workflow".
`pyfraglib` aims at being such at tool.

## Installation

`pyfraglib` can be installed using `pip`, the Python package manager:

```bash
python3 -m pip install pyfraglib
```

Since the most recent version of `pyfraglib` might not be on PyPI yet, you can
install from source, too:

```bash
git clone git@bitbucket.org:schwarzlab/project-lymphoma-cfdna.git
cd pyfraglib
python3 -m pip install .
```

## Usage

`pyfraglib` comes with a command line utility. After successful installation,
it can be used as follows:

```bash
pyfrag.py --help # show available subcommands and flags
```

`tools/workflow.sh` gives an example of a commonly used sequence of commands.

## Implemented Algorithms

| Fragmentomics Feature                    | Related Publication                 | Impl. Status  |
|------------------------------------------|-------------------------------------|---------------|
| Fragment length analysis                 |                                     | in progress   |
| K-mer end motifs 3'/5'                   |                                     | in progress   |
| Motif diversity score                    |                                     | in progress   |
| Window protection score                  |                                     | in progress   |
| D/U fragment ends ("OCF")                | Sun et al., Genome Research, 2019   | not yet       |
| Promoter fragmentation entropy           |                                     | not yet       |
| Maximum nucleosome protection            | Snyder et al., 2016                 | not yet       |
| NMF-derived fragment length signatures   | Renaud et al., Elife, 2022          | not yet       |
| DELFI Metric (short / long fragments)    | Cristiano et al., Nature 2019       | not yet       |
| 2k- & NDR-TSS coverage metrics           | Ulz et al., Nature Genetics 2016    | not yet       |

## Genomics features

Obviously, more "conventional" genomics features such as single-nucleotide
variants can be analyzed as well when dealing with cfDNA

| cfDNA Feature                         | Related Publication                     |
|---------------------------------------|-----------------------------------------|
| Genome-wide SNVs                      |                                         |
| SBS Signatures                        |                                         |
| Genome-wide CNAs                      |                                         |
| Subclonal reconstruction (ichorCNA)   | Adalsteinsson et al., Nature Comm 2017  |
| Subclonal reconstruction (liquidCNA)  | Lakatos et al., iScience 2021           |

We explicitly do not do that in our toolkit.

## License
`pyfraglib` is licensed under the GPL-v3 as indicated in the source files.
Please direct requests to <daniel.schuette@iccb-cologne.org>.

