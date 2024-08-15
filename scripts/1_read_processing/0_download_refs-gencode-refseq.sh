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
output_dir=$(realpath "${SCRIPT_DIR}/../../data/ref-gencode-refseq/")
mkdir -p "${output_dir}"

# download human reference genome and annotation
wget -N -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
wget -N -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.basic.annotation.gff3.gz"
wget -N -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/MD5SUMS"

# download pknowlesi reference genome and annotation
wget -N -P "${output_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/355/GCF_000006355.2_GCA_000006355.2/GCF_000006355.2_GCA_000006355.2_genomic.fna.gz"
wget -N -P "${output_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/355/GCF_000006355.2_GCA_000006355.2/GCF_000006355.2_GCA_000006355.2_genomic.gff.gz"
wget -N -O "${output_dir}/md5checksums-pk.txt" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/355/GCF_000006355.2_GCA_000006355.2/md5checksums.txt"

# Perform checksums
cd "${output_dir}"
sed -i '/MD5SUMS/d' "${output_dir}/MD5SUMS"
if ! md5sum --ignore-missing -c "${output_dir}/MD5SUMS" ; then
    echo "MD5 checksum failed for human reference genome..."
    exit 1
fi
if ! md5sum --ignore-missing -c "${output_dir}/md5checksums-pk.txt" ; then
    echo "MD5 checksum failed for Pk reference genome..."
    exit 1
fi

# # replace gene_biotype with gene_type attribute in Pk gff when mixing with gencode data
# (or do it the other way around since the nf-core rnaseq pipeline defaults to gene_biotype)
zcat "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.gff.gz" | sed 's/gene_biotype/gene_type/g' | gzip > "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.gene_type.gff.gz"

# concatenate human and pknowlesi genomes and annotation files
cat "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.fna.gz" \
    "${output_dir}/GRCh38.primary_assembly.genome.fa.gz" \
    >"${output_dir}/concat-GRCh38.14.gencode.46-GCF_000006355.2_GCA_000006355.2_genomic.fa.gz"

cat "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.gene_type.gff.gz" \
    "${output_dir}/gencode.v46.primary_assembly.basic.annotation.gff3.gz" \
    >"${output_dir}/concat-GRCh38.14.gencode.46-GCF_000006355.2_GCA_000006355.2_genomic.gff.gz"
