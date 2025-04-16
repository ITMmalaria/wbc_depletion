#!/bin/env bash

# set bash strict mode
set -euo pipefail

# load conda env
# conda activate wbc_depletion
export PATH="/scratch/antwerpen/203/vsc20380/containers/wbc_depletion/bin:$PATH"

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../../")

# parameter settings
star_bam_dir="${PROJECT_ROOT}/results/star_salmon"
output_dir="${PROJECT_ROOT}/results/fastq_pk_filter"
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
threads:                       ${n_threads}
"


# Extract paired fastq files from Pk-filtered bam files
for bam in "${star_bam_dir}/"*.markdup.sorted.pk.bam; do
    echo "Extracting paired fastq files from ${bam}..."
    sample_name=$(basename "${bam%_S1_L001.markdup.sorted.pk.bam}")

    samtools collate -u -O "${bam}" | \
    samtools fastq \
        -1 "${output_dir}/${sample_name}_R1.fastq.gz" \
        -2 "${output_dir}/${sample_name}_R2.fastq.gz" \
        -0 "${output_dir}/${sample_name}_R0.fastq.gz" \
        -s "${output_dir}/${sample_name}_singletons.fastq.gz" \
        -n \
        -O \
        --threads ${n_threads}
done
