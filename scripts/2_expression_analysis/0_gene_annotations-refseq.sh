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
ref_human=$(realpath "${SCRIPT_DIR}/../../data/ref-refseq-refseq/GCF_000001405.40_GRCh38.p14_genomic.gff.gz")
ref_pk=$(realpath "${SCRIPT_DIR}/../../data/ref-refseq-refseq/GCF_000006355.2_GCA_000006355.2_genomic.gff.gz")
output_file=$(realpath "${SCRIPT_DIR}/../../data/ref-refseq-refseq/gene_types.csv")

# create annotation file for rRNAs
# unlike for gencode gff files, there is no need to ignore rRNA_pseudogene lncRNA
echo "gene_type;gene_id" > "${output_file}"
zcat "${ref_human}" | awk '$3 ~ /gene/' | cut -f9 | grep -E "gene_biotype=rRNA" | grep -v "rRNA_pseudo" | grep -o "ID=[^;]*" | sed 's/ID=/human_rRNA;/' >> "${output_file}"
zcat "${ref_pk}" | cut -f9 | grep -E "gene_biotype=rRNA" | grep -o "ID=[^;]*" | sed 's/ID=/pk_rRNA;/' >> "${output_file}"
