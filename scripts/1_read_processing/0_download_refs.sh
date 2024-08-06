#!/bin/env bash

# set bash strict mode
set -euo pipefail

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

# set and create output directory
output_dir="${SCRIPT_DIR}/../../data/ref/"
mkdir -p "${output_dir}"

# download human reference genome and annotation
wget -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"
wget -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gff3.gz"
wget -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/MD5SUMS"
if [ ! md5sum --status -c "${output_dir}/MD5SUMS" ]; then
    echo "MD5 checksum failed for human reference genome."
    exit 1
fi

# download pknowlesi reference genome and annotation
wget -P "${output_dir}" https://plasmodb.org/common/downloads/release-68/PknowlesiH/fasta/data/PlasmoDB-68_PknowlesiH_Genome.fasta
wget -P "${output_dir}" "https://plasmodb.org/common/downloads/release-68/PknowlesiH/gff/data/PlasmoDB-68_PknowlesiH.gff"

# concatenate human and pknowlesi genomes
# zcat ${output_dir}/GRCh38.primary_assembly.genome.fa.gz | cat ${output_dir}/PlasmoDB-68_PknowlesiH_Genome.fasta - | bgzip >${output_dir}/concat.fa.gz
bgzip "${output_dir}/PlasmoDB-68_PknowlesiH_Genome.fasta" | cat "${output_dir}/GRCh38.primary_assembly.genome.fa.gz" >"${output_dir}/concat.fa.gz"
