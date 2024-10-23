#!/bin/bash
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
if [[ "$1" == "--clean" || "$2" == "--clean" ]]; then
    ./tools/clean.sh
    echo "Cleaned workspace."
fi

if [[ "$1" == "--install" || "$2" == "--install" ]]; then
    ./tools/local_install.sh || exit 1
    echo "Successfully installed latest \`pyfraglib' locally."
fi

SLURM_CPUS_PER_TASK=4 pyfrag.py -o frag_out/ extract --bam-dir data/ \
    --with-vcf &&
pyfrag.py -o plot_out/ stats -d frag_out/ &&
pyfrag.py -o plot_out/ lengths -d frag_out/ &&
pyfrag.py -o scores_out/ scores -d frag_out/ -b \
    data/ref/loci_covered_by_panel.bed.gz
