// This file is part of `pyfraglib`, a software suite to calculate fragmentomics
// features from cfDNA and perform downstream analyses.
//
// Copyright (C) 2024 Daniel Schütte, daniel.schuette@iccb-cologne.org
//
// This program is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version. This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details. You should have received a copy of the GNU General Public
// License along with this program. If not, see <https://www.gnu.org/licenses/>.
nextflow.enable.dsl=2

process stage_data {
    label "serial"
    tag "${sample_id}"
    errorStrategy "terminate"

    if (params.force_stage_data) { cache false }

    input:
    tuple val(sample_id), val(bam_file), val(bam_file_index), val(vcf_file)

    output:
    tuple val(sample_id), path("*.bam"), path("*.bai"), path("*.vcf")

    shell:
    """
    set +eu
    module load lang/Miniconda3
    conda activate !{params.conda_env}
    set -eu

    ln -s !{bam_file} !{sample_id}.bam
    ln -s !{bam_file_index} !{sample_id}.bam.bai
    if [[ !{vcf_file} == "none" ]]; then
        # We create a dummy VCF file to simplify downstream processing. Ugly.
        echo -e '##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' \
            > !{sample_id}.vcf
    else
        ln -s !{vcf_file} !{sample_id}.vcf
    fi
    """
}

process pyfraglib_extract {
    label "parallel"
    tag "${sample_id}"
    errorStrategy "ignore"

    if (params.force_extract) { cache false }

    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam_path), path(bai_path), path(vcf_path)

    output:
    tuple val(sample_id), path("*.frag")

    shell:
    """
    set +eu
    module load lang/Miniconda3
    conda activate !{params.conda_env}
    set -eu

    pyfrag.py -o . extract -f !{bam_path} --with-vcf
    """
}

process pyfraglib_stats {
    label "parallel"
    tag "${sample_id}"
    errorStrategy "ignore"

    if (params.force_stats) { cache false }

    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(frag_file)

    output:
    path("*.png")
    path("*.json")

    shell:
    """
    set +eu
    module load lang/Miniconda3
    conda activate !{params.conda_env}
    set -eu

    pyfrag.py -o . stats -f !{frag_file}
    """
}

process pyfraglib_lengths {
    label "parallel"
    tag "${sample_id}"
    errorStrategy { task.attempt < 4 ? "retry" : "ignore" }

    if (params.force_lengths) { cache false }

    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(frag_file)

    output:
    path("*.png")
    path("*.json")

    shell:
    """
    set +eu
    module load lang/Miniconda3
    conda activate !{params.conda_env}
    set -eu

    pyfrag.py -o . lengths -f !{frag_file} -c !{params.gmm_config_file}
    """
}

process pyfraglib_scores {
    label "huge_mem"
    tag "${sample_id}"
    errorStrategy "ignore"

    if (params.force_scores) { cache false }

    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(frag_file)

    output:
    path("*.csv")
    path("*.png")

    shell:
    """
    set +eu
    module load lang/Miniconda3
    conda activate !{params.conda_env}
    set -eu

    pyfrag.py -o . scores -f !{frag_file} -b !{params.bed_file}
    """
}

workflow {
    sample_info = Channel.fromPath(params.sample_tbl)
        | splitCsv(header: false, sep: '\t')
        | map { row -> tuple(row[0], row[1], row[2], row[3]) }
    samples = stage_data(sample_info)

    frag_files = pyfraglib_extract(samples)

    pyfraglib_stats(frag_files)
    pyfraglib_lengths(frag_files)
    pyfraglib_scores(frag_files)
}
