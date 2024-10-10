#!/bin/env bash

# set bash strict mode
set -euo pipefail

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../../")

# parameter settings
samplesheet="${PROJECT_ROOT}/data/samplesheet.csv"
output="${PROJECT_ROOT}/results-refseq/mapping_stats.csv"

# echo "sample,condition,RNA_type,library_kit,removal,WBC_depletion,kit,primary_mapped_all,primary_mapped_pk,fastq_screen_pk_onehit,fastq_screen_pk_multihit" > "${output}"
echo "sample,condition,RNA_type,library_kit,removal,WBC_depletion,kit,primary_mapped_all,primary_mapped_pk" > "${output}"

tail -n +2 "${samplesheet}" | \
    while IFS=, read -r sample fastq_1 fastq_2 strandedness condition RNA_type library_kit removal WBC_depletion kit; do
        # extract properly paired count from samtools flagstat
        primary_mapped_all=$(sed '8q;d' "${PROJECT_ROOT}/results-refseq/star_salmon/samtools_stats/${sample}.markdup.sorted.bam.flagstat" | sed 's/ \+.*//g')
        primary_mapped_pk=$(sed '8q;d' "${PROJECT_ROOT}/results-refseq/star_salmon/samtools_stats_pk/${sample}.markdup.sorted.pk.bam.flagstat" | sed 's/ \+.*//g')

#         # extract properly paired count from samtools flagstat
#         proper_mapped_all=$(sed '12q;d' "${PROJECT_ROOT}/results-refseq/star_salmon/samtools_stats/${sample}.markdup.sorted.bam.flagstat" | sed 's/ \+.*//g')
#         proper_mapped_pk=$(sed '12q;d' "${PROJECT_ROOT}/results-refseq/star_salmon/samtools_stats_pk/${sample}.markdup.sorted.pk.bam.flagstat" | sed 's/ \+.*//g')

#         # extract fastq screen info - now read directly in R
#         fastq_screen_pk_onehit_r1=$(grep "PknowlesiH" "${PROJECT_ROOT}/results-refseq/fastq-screen/${sample}_R1_001_screen.txt" | cut -f5)
#         fastq_screen_pk_onehit_r2=$(grep "PknowlesiH" "${PROJECT_ROOT}/results-refseq/fastq-screen/${sample}_R2_001_screen.txt" | cut -f5)
#         fastq_screen_pk_multihit_r1=$(grep "PknowlesiH" "${PROJECT_ROOT}/results-refseq/fastq-screen/${sample}_R1_001_screen.txt" | cut -f7)
#         fastq_screen_pk_multihit_r2=$(grep "PknowlesiH" "${PROJECT_ROOT}/results-refseq/fastq-screen/${sample}_R2_001_screen.txt" | cut -f7)

#         echo ${sample},${condition},${RNA_type},${library_kit},${removal},${WBC_depletion},${kit},${primary_mapped_all},${primary_mapped_pk},$((${fastq_screen_pk_onehit_r1} + ${fastq_screen_pk_onehit_r2})),$((${fastq_screen_pk_multihit_r1}+${fastq_screen_pk_multihit_r2})) >> "${output}"
        echo ${sample},${condition},${RNA_type},${library_kit},${removal},${WBC_depletion},${kit},${primary_mapped_all},${primary_mapped_pk} >> "${output}"
#         echo ${sample},${condition},${RNA_type},${library_kit},${removal},${WBC_depletion},${kit},${proper_mapped_all},${proper_mapped_pk} >> "${output}"
    done
