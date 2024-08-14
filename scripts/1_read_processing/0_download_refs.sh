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
# wget -N -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"
# wget -N -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gff3.gz"
# wget -N -P "${output_dir}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/MD5SUMS"
wget -N -P "${output_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
wget -N -P "${output_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
wget -N -P "${output_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/md5checksums.txt"

# download pknowlesi reference genome and annotation
# wget -nc -P "${output_dir}" "https://plasmodb.org/common/downloads/release-68/PknowlesiH/fasta/data/PlasmoDB-68_PknowlesiH_Genome.fasta"
# wget -nc -P "${output_dir}" "https://plasmodb.org/common/downloads/release-68/PknowlesiH/gff/data/PlasmoDB-68_PknowlesiH.gff"
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
    >"${output_dir}/concat_GRCh38.14_PkH.fa.gz"

cat "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.gff.gz" \
    "${output_dir}/GCF_000001405.40_GRCh38.p14_genomic.gff.gz" \
    >"${output_dir}/concat_GRCh38.14_PkH.gff.gz"

############

# gencode human - refseq pk

# # replace gene_biotype with gene_type attribute in Pk gff when mixing with gencode data
# (or do it the other way around since the nf-core rnaseq pipeline defaults to gene_biotype)
# zcat "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.gff.gz" | sed 's/gene_biotype/gene_type/g' | gzip > "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.gene_type.gff.gz"

# # concatenate human and pknowlesi genomes and annotation files
# cat "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.fna.gz" "${output_dir}/GRCh38.primary_assembly.genome.fa.gz" >"${output_dir}/concat_GRCh38-v45_GCF_000006355.2.fa.gz"
# cat "${output_dir}/GCF_000006355.2_GCA_000006355.2_genomic.gene_type.gff.gz" "${output_dir}/gencode.v45.annotation.gff3.gz" >"${output_dir}/concat_GRCh38-v45_GCF_000006355.2.gff.gz"

############

# gencode human - plasmodb pk

# activate conda env with samtools installed
# conda create -n rnaseq
# conda install -c bioconda samtools
# source /data/antwerpen/203/vsc20380/miniforge3/etc/profile.d/conda.sh
# conda activate rnaseq
# uses bgzip as bundled with samtools, see https://www.htslib.org/doc/bgzip.html
## do not use pipes to ensure bgzip conversion is succesful, otherwise empty concat file is created
## bgzip --keep "${output_dir}/PlasmoDB-68_PknowlesiH_Genome.fasta" | cat "${output_dir}/GRCh38.primary_assembly.genome.fa.gz" - >"${output_dir}/concat_GRCh38-v45_PkH-PlasmoDB-v68.fa.gz"
## bgzip --keep "${output_dir}/PlasmoDB-68_PknowlesiH.gff" | cat "${output_dir}/gencode.v45.annotation.gff3.gz" - >"${output_dir}/concat_GRCh38-v45_PkH-PlasmoDB-v68.gff.gz"
# bgzip --keep --force "${output_dir}/PlasmoDB-68_PknowlesiH_Genome.fasta"
# bgzip --keep --force "${output_dir}/PlasmoDB-68_PknowlesiH.gff"
# cat "${output_dir}/PlasmoDB-68_PknowlesiH_Genome.fasta.gz" "${output_dir}/GRCh38.primary_assembly.genome.fa.gz" >"${output_dir}/concat_GRCh38-v45_PkH-PlasmoDB-v68.fa.gz"
# cat "${output_dir}/PlasmoDB-68_PknowlesiH.gff.gz" "${output_dir}/gencode.v45.annotation.gff3.gz" >"${output_dir}/concat_GRCh38-v45_PkH-PlasmoDB-v68.gff.gz"
