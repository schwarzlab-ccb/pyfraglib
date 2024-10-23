#!/bin/sh
# This pipeline is less than sophisticated and was created within an hour.
# Do not use for real tasks!
#
# `bcftools' and `snpEff' need to be installed and an existing BAM file must
# be provided as an input. A VCF file of the same name will be created. We do
# only very limited filtering to not loose too many variants. That's reasonable
# for testing, but certainly not for real science.
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
if [ ! -e "$1" ]; then
    printf "No input file provided or \`%s' does not exist.\n" "$1"
    exit 1
fi
infile="$1"
basename="${infile%.*}"

printf "Annotating variants in \`%s'.\n" "$infile"
bcftools mpileup -f ref/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    "$infile" -Ou -o temp_out.vcf
bcftools call -mv -Ov -o temp_out1.vcf temp_out.vcf
bcftools +fill-tags temp_out1.vcf -- -t AF > temp_out2.vcf
bcftools view -v snps temp_out2.vcf -o temp_out3.vcf

# Allele frequencies of true tumor variants tend to be very low.
bcftools filter -i 'QUAL>30 & DP>=1 & DP<2000 & AF>0.0001 & AF<0.95' \
    temp_out3.vcf -o temp_out4.vcf
snpEff GRCh37.p13 temp_out4.vcf > "$basename".vcf

# Clean up build artifacts.
echo "Cleaning up."
rm -f snpEff_genes.txt snpEff_summary.html temp_out*
