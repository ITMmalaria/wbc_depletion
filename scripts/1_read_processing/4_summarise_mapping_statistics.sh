#!/bin/env bash

# set bash strict mode
set -euo pipefail

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/../../")

# parameter settings
samplesheet="${PROJECT_ROOT}/data/samplesheet.csv"
output="${PROJECT_ROOT}/results/mapping_stats.csv"

echo "sample,condition,RNA_type,library_kit,removal,WBC_depletion,kit,trimgalore_total_1,trimgalore_filtered_1,trimgalore_total_2,trimgalore_filtered_2,star_input_reads,star_unique_mapped_n,star_unique_mapped_pct,star_splices_n,star_multiple_loci_n,star_multiple_loci_pct,star_too_many_loci_n,star_too_many_loci_pct,star_unmapped_mismatch_n,star_unmapped_mismatch_pct,star_unmapped_short_n,star_unmapped_short_pct,star_unmapped_other_n,star_unmapped_other_pct,star_chimeric_n,star_chimeric_pct,unique_mapped_pk,flagstat_total_all,flagstat_total_pk,flagstat_primary_all,flagstat_primary_pk,flagstat_secondary_all,flagstat_secondary_pk,flagstat_primary_mapped_all,flagstat_primary_mapped_pk,flagstat_proper_pair,flagstat_proper_pair_pk" > "${output}"

tail -n +2 "${samplesheet}" | \
    while IFS=, read -r sample fastq_1 fastq_2 strandedness condition RNA_type library_kit removal WBC_depletion kit; do

        # extract stats from trimgalore
        trim_log="${PROJECT_ROOT}/results/trimgalore/${sample}_1.fastq.gz_trimming_report.txt"
        trimgalore_total_1=$(grep "Total reads processed" ${trim_log} | awk '{ print $NF }' | sed 's/,//g')
        trimgalore_filtered_1=$(grep "Reads written (passing filters)" ${trim_log} | awk '{ print $(NF-1) }' | sed 's/,//g')
        trim_log="${PROJECT_ROOT}/results/trimgalore/${sample}_2.fastq.gz_trimming_report.txt"
        trimgalore_total_2=$(grep "Total reads processed" ${trim_log} | awk '{ print $NF }' | sed 's/,//g')
        trimgalore_filtered_2=$(grep "Reads written (passing filters)" ${trim_log} | awk '{ print $(NF-1) }' | sed 's/,//g')

        # extract stats from STAR .Log.final.out
        log="${PROJECT_ROOT}/results/star_salmon/log/${sample}.Log.final.out"

        star_input_reads=$(grep "Number of input reads" ${log} | cut -d '|' -f2 | tr -d '[:space:]')

        star_unique_mapped_n=$(grep "Uniquely mapped reads number" ${log} | cut -d '|' -f2 | tr -d '[:space:]')
        star_unique_mapped_pct=$(grep "Uniquely mapped reads %" ${log} | cut -d '|' -f2 | tr -d '[:space:]')

        star_splices_n=$(grep "Number of splices: Total" ${log} | cut -d '|' -f2 | tr -d '[:space:]')

        star_multiple_loci_n=$(grep "Number of reads mapped to multiple loci" ${log} | cut -d '|' -f2 | tr -d '[:space:]')
        star_multiple_loci_pct=$(grep "% of reads mapped to multiple loci" ${log} | cut -d '|' -f2 | tr -d '[:space:]')

        star_too_many_loci_n=$(grep "Number of reads mapped to too many loci " ${log} | cut -d '|' -f2 | tr -d '[:space:]')
        star_too_many_loci_pct=$(grep "% of reads mapped to too many loci" ${log} | cut -d '|' -f2 | tr -d '[:space:]')

        star_unmapped_mismatch_n=$(grep "Number of reads unmapped: too many mismatches" ${log} | cut -d '|' -f2 | tr -d '[:space:]')
        star_unmapped_mismatch_pct=$(grep "% of reads unmapped: too many mismatches" ${log} | cut -d '|' -f2 | tr -d '[:space:]')

        star_unmapped_short_n=$(grep "Number of reads unmapped: too short" ${log} | cut -d '|' -f2 | tr -d '[:space:]')
        star_unmapped_short_pct=$(grep "% of reads unmapped: too short" ${log} | cut -d '|' -f2 | tr -d '[:space:]')

        star_unmapped_other_n=$(grep "Number of reads unmapped: other" ${log} | cut -d '|' -f2 | tr -d '[:space:]')
        star_unmapped_other_pct=$(grep "% of reads unmapped: other" ${log} | cut -d '|' -f2 | tr -d '[:space:]')

        star_chimeric_n=$(grep "Number of chimeric reads" ${log} | cut -d '|' -f2 | tr -d '[:space:]')
        star_chimeric_pct=$(grep "% of chimeric reads" ${log} | cut -d '|' -f2 | tr -d '[:space:]')

        # count uniquely mapped pk reads (same calculation as STAR uniquely mapped)
        # see: https://github.com/alexdobin/STAR/issues/507
        unique_mapped_pk=$(samtools view "${PROJECT_ROOT}/results/star_salmon/${sample}.markdup.sorted.pk.bam" | awk '/NH:i:1\t/ {print $1} ' | sort | uniq | wc -l)

        # extract properly paired count from samtools flagstat
        flagstat_total_all=$(sed '1q;d' "${PROJECT_ROOT}/results/star_salmon/samtools_stats/${sample}.markdup.sorted.bam.flagstat" | sed 's/ \+.*//g')
        flagstat_total_pk=$(sed '1q;d' "${PROJECT_ROOT}/results/star_salmon/samtools_stats_pk/${sample}.markdup.sorted.pk.bam.flagstat" | sed 's/ \+.*//g')

        flagstat_primary_all=$(sed '2q;d' "${PROJECT_ROOT}/results/star_salmon/samtools_stats/${sample}.markdup.sorted.bam.flagstat" | sed 's/ \+.*//g')
        flagstat_primary_pk=$(sed '2q;d' "${PROJECT_ROOT}/results/star_salmon/samtools_stats_pk/${sample}.markdup.sorted.pk.bam.flagstat" | sed 's/ \+.*//g')

        flagstat_secondary_all=$(sed '3q;d' "${PROJECT_ROOT}/results/star_salmon/samtools_stats/${sample}.markdup.sorted.bam.flagstat" | sed 's/ \+.*//g')
        flagstat_secondary_pk=$(sed '3q;d' "${PROJECT_ROOT}/results/star_salmon/samtools_stats_pk/${sample}.markdup.sorted.pk.bam.flagstat" | sed 's/ \+.*//g')

        flagstat_primary_mapped_all=$(sed '8q;d' "${PROJECT_ROOT}/results/star_salmon/samtools_stats/${sample}.markdup.sorted.bam.flagstat" | sed 's/ \+.*//g')
        flagstat_primary_mapped_pk=$(sed '8q;d' "${PROJECT_ROOT}/results/star_salmon/samtools_stats_pk/${sample}.markdup.sorted.pk.bam.flagstat" | sed 's/ \+.*//g')

        flagstat_proper_pair=$(sed '12q;d' "${PROJECT_ROOT}/results/star_salmon/samtools_stats/${sample}.markdup.sorted.bam.flagstat" | sed 's/ \+.*//g')
        flagstat_proper_pair_pk=$(sed '12q;d' "${PROJECT_ROOT}/results/star_salmon/samtools_stats_pk/${sample}.markdup.sorted.pk.bam.flagstat" | sed 's/ \+.*//g')

    echo ${sample},${condition},${RNA_type},${library_kit},${removal},${WBC_depletion},${kit},${trimgalore_total_1},${trimgalore_filtered_1},${trimgalore_total_2},${trimgalore_filtered_2},${star_input_reads},${star_unique_mapped_n},${star_unique_mapped_pct},${star_splices_n},${star_multiple_loci_n},${star_multiple_loci_pct},${star_too_many_loci_n},${star_too_many_loci_pct},${star_unmapped_mismatch_n},${star_unmapped_mismatch_pct},${star_unmapped_short_n},${star_unmapped_short_pct},${star_unmapped_other_n},${star_unmapped_other_pct},${star_chimeric_n},${star_chimeric_pct},${unique_mapped_pk},${flagstat_total_all},${flagstat_total_pk},${flagstat_primary_all},${flagstat_primary_pk},${flagstat_secondary_all},${flagstat_secondary_pk},${flagstat_primary_mapped_all},${flagstat_primary_mapped_pk},${flagstat_proper_pair},${flagstat_proper_pair_pk} >> "${output}"

    done
