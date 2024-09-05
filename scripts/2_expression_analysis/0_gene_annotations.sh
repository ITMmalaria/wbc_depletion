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

# set input and output
ref_human=$(realpath "${SCRIPT_DIR}/../../data/ref/gencode.v46.primary_assembly.basic.annotation.gff3.gz")
ref_pk=$(realpath "${SCRIPT_DIR}/../../data/ref/PlasmoDB-68_PknowlesiH.gff")
output_file=$(realpath "${SCRIPT_DIR}/../../data/ref/gene_types.csv")

# create annotation file for rRNAs
echo "gene_type;gene_id" > "${output_file}"
zcat "${ref_human}"| awk '$3 ~ /gene/' | cut -f9 | grep -E "gene_type=rRNA" | grep -v "rRNA_pseudo" | grep -o "ID=[^;]*" | sed 's/ID=/human_rRNA;/' >> "${output_file}"
zcat "${ref_human}"| awk '$3 ~ /gene/' | cut -f9 | grep -E "gene_type=Mt_rRNA" | grep -o "ID=[^;]*" | sed 's/ID=/human_Mt_rRNA;/' >> "${output_file}"
cut -f9 "${ref_pk}" | grep ";ebi_biotype=rRNA" | grep -o "ID=[^;]*" | sed 's/ID=/pk_rRNA;/' >> "${output_file}"

# note: gene_ebi_biotype=rRNA entries are categorised as rRNA
# grep "gene_ebi_biotype=rRNA" data/ref/PlasmoDB-68_PknowlesiH.gff | cut -f3 | sort -u
# rRNA
# while ebi_biotype=rRNA entries are ncRNA_genes
# grep ";ebi_biotype=rRNA" data/ref/PlasmoDB-68_PknowlesiH.gff | cut -f3 | sort -u
# ncRNA_gene
# but the former are all transcripts with the latter as their parent
