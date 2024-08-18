#!/bin/env bash

# set bash strict mode
set -euo pipefail

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')")
else
    # SCRIPT_DIR=$(realpath "$0")
    SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
fi

# set and create output directory
output_dir=$(realpath "${SCRIPT_DIR}/../../data/ref/")
mkdir -p "${output_dir}"

# download human reference genome and annotation
wget -N -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
wget -N -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.basic.annotation.gff3.gz"
wget -N -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/MD5SUMS"

# download pknowlesi reference genome and annotation
wget -nc -P "${output_dir}" "https://plasmodb.org/common/downloads/release-68/PknowlesiH/fasta/data/PlasmoDB-68_PknowlesiH_Genome.fasta"
wget -nc -P "${output_dir}" "https://plasmodb.org/common/downloads/release-68/PknowlesiH/gff/data/PlasmoDB-68_PknowlesiH.gff"

# Perform checksums
cd "${output_dir}"
sed -i '/MD5SUMS/d' "${output_dir}/MD5SUMS"
if ! md5sum --ignore-missing -c "${output_dir}/MD5SUMS"; then
    echo "MD5 checksum failed for human reference genome..."
    exit 1
fi

# concatenate human and pknowlesi genomes and annotation files
gzip -c "${output_dir}/PlasmoDB-68_PknowlesiH_Genome.fasta" |
    cat "${output_dir}/GRCh38.primary_assembly.genome.fa.gz" - \
        >"${output_dir}/concat_GRCh38.14.gencode.46_PkH.PlasmoDB.68.fa.gz"

gzip -c "${output_dir}/PlasmoDB-68_PknowlesiH.gff" |
    cat "${output_dir}/gencode.v46.primary_assembly.basic.annotation.gff3.gz" - \
        >"${output_dir}/concat_GRCh38.14.gencode.46_PkH.PlasmoDB.68.gff.gz"
