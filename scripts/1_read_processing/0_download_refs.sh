#!/bin/env bash

# set bash strict mode
set -euo pipefail

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
if [ -n "${SLURM_JOB_ID:-}" ] ; then
    SCRIPT_DIR=$(dirname $(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}'))
else
    # SCRIPT_DIR=$(realpath "$0")
    SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
fi

# set and create output directory
output_dir=$(realpath "${SCRIPT_DIR}/../../data/ref/")
mkdir -p "${output_dir}"

# download human reference genome and annotation
wget -N -P "${output_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2023_10/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
wget -N -P "${output_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2023_10/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
wget -N -P "${output_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2023_10/md5checksums.txt"

# download pknowlesi reference genome and annotation
wget -N -P "${output_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/355/GCF_000006355.2_GCA_000006355.2/GCF_000006355.2_GCA_000006355.2_genomic.fna.gz"
wget -N -P "${output_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/355/GCF_000006355.2_GCA_000006355.2/GCF_000006355.2_GCA_000006355.2_genomic.gff.gz"
wget -N -O "${output_dir}/md5checksums-pk.txt" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/355/GCF_000006355.2_GCA_000006355.2/md5checksums.txt"

# Perform checksums
cd "${output_dir}"
if ! md5sum --ignore-missing -c "${output_dir}/md5checksums.txt" ; then
    echo "MD5 checksum failed for human reference genome..."
    exit 1
fi
if ! md5sum --ignore-missing -c "${output_dir}/md5checksums-pk.txt" ; then
    echo "MD5 checksum failed for Pk reference genome..."
    exit 1
fi

# concatenate human and pknowlesi genomes and annotation files
cat "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.fna.gz" \
    "${output_dir}/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" \
    >"${output_dir}/concat-GCF_000001405.40_GRCh38.p14_genomic-GCF_000006355.2_GCA_000006355.2_genomic.fa.gz"

cat "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.gff.gz" \
    "${output_dir}/GCF_000001405.40_GRCh38.p14_genomic.gff.gz" \
    >"${output_dir}/concat-GCF_000001405.40_GRCh38.p14_genomic-GCF_000006355.2_GCA_000006355.2_genomic.gff.gz"
