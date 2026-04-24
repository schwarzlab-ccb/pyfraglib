// This file is part of `pyfraglib`, a software suite to calculate fragmentomics
// features from cfDNA and perform downstream analyses.
//
// Copyright (C) 2026 Daniel Schütte, daniel.schuette@iccb-cologne.org
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
    tuple val(sample_id), path("*.frag"), emit: frag
    path("*_extract_stats.tsv"), emit: stats

    shell:
    """
    set +eu
    module load lang/Miniconda3
    conda activate !{params.conda_env}
    set -eu

    bam_bytes=\$(stat -L -c %s !{bam_path})
    bam_reads_total=\$(samtools view -c !{bam_path})
    bam_reads_dedup=\$(samtools view -c -F 1024 !{bam_path})

    /usr/bin/time -f '%e\\t%U\\t%S\\t%M' -o extract.time \\
        pyfrag.py extract -o . -f !{bam_path} --with-vcf

    extract_wall_s=\$(awk -F'\\t' '{print \$1}' extract.time)
    extract_user_s=\$(awk -F'\\t' '{print \$2}' extract.time)
    extract_sys_s=\$(awk -F'\\t' '{print \$3}' extract.time)
    extract_max_rss_kb=\$(awk -F'\\t' '{print \$4}' extract.time)

    frag_path=\$(ls *.frag | head -n 1)
    frag_bytes=\$(stat -c %s "\$frag_path")
    frag_fragments=\$(python -c "
import sys
import pyarrow.parquet as pq
print(pq.ParquetFile(sys.argv[1]).metadata.num_rows)
" "\$frag_path")

    compression_ratio_total=\$(awk -v b=\$bam_bytes -v f=\$frag_bytes 'BEGIN {if (f > 0) printf "%.4f", b / f; else print "nan"}')
    compression_ratio_dedup=\$(awk -v b=\$bam_bytes -v rt=\$bam_reads_total -v rd=\$bam_reads_dedup -v f=\$frag_bytes 'BEGIN {if (f > 0 && rt > 0) printf "%.4f", (b * rd / rt) / f; else print "nan"}')
    reads_per_sec=\$(awk -v n=\$bam_reads_total -v t=\$extract_wall_s 'BEGIN {if (t > 0) printf "%.1f", n / t; else print "nan"}')
    fragments_per_sec=\$(awk -v n=\$frag_fragments -v t=\$extract_wall_s 'BEGIN {if (t > 0) printf "%.1f", n / t; else print "nan"}')

    stats_file="!{sample_id}_extract_stats.tsv"
    printf 'sample_id\\tbam_bytes\\tbam_reads_total\\tbam_reads_dedup\\tfrag_bytes\\tfrag_fragments\\textract_wall_s\\textract_user_s\\textract_sys_s\\textract_max_rss_kb\\treads_per_sec\\tfragments_per_sec\\tcompression_ratio_total\\tcompression_ratio_dedup\\n' > "\$stats_file"
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "!{sample_id}" "\$bam_bytes" "\$bam_reads_total" "\$bam_reads_dedup" \\
        "\$frag_bytes" "\$frag_fragments" \\
        "\$extract_wall_s" "\$extract_user_s" "\$extract_sys_s" "\$extract_max_rss_kb" \\
        "\$reads_per_sec" "\$fragments_per_sec" \\
        "\$compression_ratio_total" "\$compression_ratio_dedup" >> "\$stats_file"
    """
}


process collect_extract_stats {
    label "serial"
    tag "cohort"
    errorStrategy "terminate"

    publishDir "${params.out_dir}", mode: 'copy', overwrite: true

    input:
    path(per_sample_stats)

    output:
    path("extract_stats.tsv")

    shell:
    """
    first=\$(ls !{per_sample_stats} | head -n 1)
    head -n 1 "\$first" > extract_stats.tsv
    for f in !{per_sample_stats}; do
        tail -n +2 "\$f"
    done | sort -u >> extract_stats.tsv
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
    path("*.csv")

    shell:
    """
    set +eu
    module load lang/Miniconda3
    conda activate !{params.conda_env}
    set -eu

    pyfrag.py stats -o . -f !{frag_file}
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

    pyfrag.py lengths -o . -f !{frag_file} -c !{params.gmm_config_file}
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

    pyfrag.py scores -o . -f !{frag_file} -b !{params.bed_file}
    """
}

workflow {
    sample_info = Channel.fromPath(params.sample_tbl)
        | splitCsv(header: false, sep: '\t')
        | map { row -> tuple(row[0], row[1], row[2], row[3]) }
    samples = stage_data(sample_info)

    extract_result = pyfraglib_extract(samples)
    frag_files = extract_result.frag

    pyfraglib_stats(frag_files)
    pyfraglib_lengths(frag_files)
    pyfraglib_scores(frag_files)
    collect_extract_stats(extract_result.stats.collect())
}
