#!/bin/bash -l

#SBATCH --job-name=gencode-refseq-nfcore-rnaseq      # create a short name for your job
#SBATCH --nodes=1                                   # node count
#SBATCH --ntasks=1                                  # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=12:00:00                             # total run time limit (HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL                  # send email when job begins, ends or fails
#SBATCH --mail-user=pmoris@itg.be
#SBATCH -o slurm-stdout.%j.%x.out
#SBATCH -e slurm-stderr.%j.%x.out

# load Java module on HPC to enable Nextflow
module load Java

printf "Script execution started at $(date).\n"

# get file path of script
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
# if [ -n "${SLURM_JOB_ID:-}" ] ; then
    # PROJECT_ROOT=$(dirname $(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}'))
# else
    # PROJECT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
# fi
PROJECT_ROOT="/scratch/antwerpen/203/vsc20380/wbc_depletion"

nextflow run nf-core/rnaseq \
    -r 3.14.0 \
    --input "${PROJECT_ROOT}/data/samplesheet-subset.csv" \
    --outdir "${PROJECT_ROOT}/results-gencode-refseq-onepass" \
    --fasta "${PROJECT_ROOT}/data/ref-gencode-refseq/concat-GRCh38.14.gencode.46-GCF_000006355.2_GCA_000006355.2_genomic.fa.gz" \
    --gff "${PROJECT_ROOT}/data/ref-gencode-refseq/concat-GRCh38.14.gencode.46-GCF_000006355.2_GCA_000006355.2_genomic.gff.gz" \
    --aligner star_salmon \
    --trimmer trimgalore \
    --extra_salmon_quant_args '--seqBias --gcBias' \
    -c "${PROJECT_ROOT}/scripts/1_read_processing/STAR_ALIGN-onepass.config" \
    -profile vsc_calcua \
    -ansi-log false \
    --skip_biotype_qc \
    --save_reference \
    -resume

    # --skip_alignment \
    # --skip_pseudo_alignment \
