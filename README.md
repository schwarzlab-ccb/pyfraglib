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

It is recommended that `pyfraglib` is installed into a dedicated conda
environment using `pip`, the Python package manager:

```bash
git clone git@bitbucket.org:schwarzlab/project-lymphoma-cfdna.git
cd pyfraglib
conda env create -f pyfraglib.yml
conda activate pyfraglib
python3 -m pip install .
```

During development, all code in `pyfraglib` is thoroughly type-checked using
`mypy`. After an initial installation as described above, do:

```bash
./tools/dev_install.sh # typing & linting errors are reported
```

As soon as `pyfraglib` is available through PyPI, installation will be even
easier (no need to clone the repo then).

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
| NMF-derived fragment length signatures   | Renaud et al., Elife, 2022          | not yet       |
| DELFI Metric (short / long fragments)    | Cristiano et al., Nature 2019       | not yet       |
| Maximum nucleosome protection            | Snyder et al., 2016                 | not planned   |
| 2k- & NDR-TSS coverage metrics           | Ulz et al., Nature Genetics 2016    | not planned   |
| Promoter fragmentation entropy           |                                     | not planned   |

## License
`pyfraglib` is licensed under the GPL-v3 as indicated in the source files.
Please direct requests to <daniel.schuette@iccb-cologne.org>.

