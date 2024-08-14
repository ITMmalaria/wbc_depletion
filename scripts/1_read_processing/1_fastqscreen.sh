#!/bin/env bash

# set bash strict mode
set -euo pipefail

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
PROJECT_DIR=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../../")

# parameter settings
fastq_dir="${PROJECT_DIR}/data/fastq"
output_dir="${PROJECT_DIR}/results/fastq-screen"
fastq_screen_conf="${PROJECT_DIR}/config/fastq-screen.conf"
ref_human="${PROJECT_DIR}/data/ref-gencode-plasmodb/GRCh38.primary_assembly.genome.fa.gz"
ref_pk="${PROJECT_DIR}/data/ref-gencode-plasmodb/PlasmoDB-68_PknowlesiH_Genome.fasta"
n_threads="${SLURM_CPUS_PER_TASK:-8}"

# create output directories
mkdir -p "${output_dir}"

# check if reference fasta files exists
for ref in ${ref_human} ${ref_pk}; do
    if ! [ -f "${ref}" ]; then
        echo "Reference fasta file not found (${ref})."
        exit 1
    fi
done

# log run options
printf "
FastQ Screen script | $(basename "$0")
==============================================

Output directory:           ${output_dir}
FASTQ reads directory:      ${fastq_dir}
Reference human:            ${ref_human}
Reference P. knowlesi:      ${ref_pk}
FastQ Screen config file:   ${fastq_screen_conf}
threads:                    ${n_threads}
"

# create reference index if it does not yet exist
# required for fastq-screen
for ref in ${ref_human} ${ref_pk}; do
    for i in "${ref}."{amb,ann,bwt,pac,sa}; do
        if ! [ -f "${i}" ]; then
            index_files_found=0
            echo "Building BWA index for ${ref}..."
            bwa index "${ref}"
            break
        else
            index_files_found=1
        fi
    done
    if [ "$index_files_found" -eq 1 ]; then
        echo "Found BWA index files for ${ref}, skipping indexing step..."
    fi
done

# run FastQ Screen (threads option is inherited by bwa/bowtie)
for read in "${fastq_dir}/"*.fastq.gz; do
    fastq_screen \
        --threads "${n_threads}" \
        --aligner bwa \
        --conf "${fastq_screen_conf}" \
        --outdir "${output_dir}" \
        "${read}"
done
