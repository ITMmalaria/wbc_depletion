#!/bin/env bash

# set bash strict mode
set -euo pipefail

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../../")

# parameter settings
star_bam_dir="${PROJECT_ROOT}/results/star_salmon"
output_dir="${PROJECT_ROOT}/results/star_salmon/samtools_stats_pk"
ref_pk="${PROJECT_ROOT}/data/ref-refseq-refseq/GCF_000006355.2_GCA_000006355.2_genomic.fna"
gff_pk="${PROJECT_ROOT}/data/ref-refseq-refseq/GCF_000006355.2_GCA_000006355.2_genomic.gff.gz"

n_threads=$(("${SLURM_CPUS_PER_TASK:-8}" - 1))

# create output directories
mkdir -p "${output_dir}"

echo ${output_dir}

# log run options
printf "
Samtools filter and flagstat script script | $(basename "$0")
==============================================

Output directory:              ${output_dir}
STAR .bam directory:           ${star_bam_dir}
P. knowlesi reference fasta:   ${ref_pk}
P. knowlesi annotation gff:    ${gff_pk}
threads:                       ${n_threads}
"

# Create bed file with chromsomes/regions based on Pk annotation gff file.
samtools faidx "${ref_pk}"
bed_pk="${ref_pk%.fasta}.bed"
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' "${ref_pk}.fai" >"${bed_pk}"

# Filter on Pk
for bam in "${star_bam_dir}/"*.markdup.sorted.bam; do
    echo "Processing file ${bam}..."
    pk_bam=$(basename "${bam%.bam}.pk.bam")
    samtools view --threads "${n_threads}" \
        -b -h -L "${bed_pk}" \
        "${bam}" \
        >"${star_bam_dir}/${pk_bam}"
    samtools index --threads "${n_threads}" \
        "${star_bam_dir}/${pk_bam}"
    samtools stats --threads "${n_threads}" \
        "${star_bam_dir}/${pk_bam}" \
        >"${output_dir}/${pk_bam}.stats"
    samtools flagstat --threads "${n_threads}" \
        "${star_bam_dir}/${pk_bam}" \
        >"${output_dir}/${pk_bam}.flagstat"
    samtools idxstats --threads "${n_threads}" \
        "${star_bam_dir}/${pk_bam}" \
        >"${output_dir}/${pk_bam}.idxstats"
done

# # Filter on human
# for bam in "${star_bam_dir}/"*.markdup.sorted.bam; do
#     echo "Processing file ${bam}..."
#     human_bam=$(basename "${bam%.bam}.human.bam")
#     samtools view --threads "${n_threads}" -b -h -L -U "${bed_pk}" "${bam}" >"${star_bam_dir}/${human_bam}"
#     samtools index --threads "${n_threads}" "${human_bam}"
#     samtools stats --threads "${n_threads}" "${human_bam}" >"${output_dir}/${human_bam}.stats"
#     samtools flagstat --threads "${n_threads}" "${human_bam}" >"${output_dir}/${human_bam}.flagstat"
#     samtools idxstats --threads "${n_threads}" "${human_bam}" >"${output_dir}/${human_bam}.idxstats"
# done
