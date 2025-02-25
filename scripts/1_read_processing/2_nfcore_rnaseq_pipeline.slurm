#!/bin/bash -l

#SBATCH --job-name=nfcore-rnaseq                    # create a short name for your job
#SBATCH --nodes=1                                   # node count
#SBATCH --ntasks=1                                  # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=24:00:00                             # total run time limit (HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL                  # send email when job begins, ends or fails
#SBATCH --mail-user=pmoris@itg.be
#SBATCH -o slurm-stdout.%j.%x.out
#SBATCH -e slurm-stderr.%j.%x.out

# load Java module on HPC to enable Nextflow (alternatively load conda environment with Nextflow)
module load Java

printf "Script execution started at $(date).\n"

# Get file path of script and set project root.
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')")
else
    SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
fi
PROJECT_ROOT=$(realpath "${SCRIPT_DIR}/../../")

# version 3.14.0 of the nf-core/rnaseq pipeline does not yet allow fixed parameters
# to be overridden via the cli, therefor a custom config file is required instead
# see: https://github.com/nf-core/rnaseq/issues/1046
# --extra_star_align_args '--twopassMode None' \
nextflow run nf-core/rnaseq \
-r 3.14.0 \
--input "${PROJECT_ROOT}/data/samplesheet.csv" \
--outdir "${PROJECT_ROOT}/results" \
--fasta "${PROJECT_ROOT}/data/ref/concat-GCF_000001405.40_GRCh38.p14_genomic-GCF_000006355.2_GCA_000006355.2_genomic.fa.gz" \
--gff "${PROJECT_ROOT}/data/ref/concat-GCF_000001405.40_GRCh38.p14_genomic-GCF_000006355.2_GCA_000006355.2_genomic.gff.gz" \
--aligner star_salmon \
--trimmer trimgalore \
--extra_salmon_quant_args '--seqBias --gcBias' \
--skip_biotype_qc \
-c "${PROJECT_ROOT}/config/STAR_ALIGN-onepass.config" \
-profile vsc_calcua \
-ansi-log false \
--save_reference \
-resume

# for resuming or re-running, --save-reference can be omitted and replaced with
# --star_index "${PROJECT_ROOT}/results/genome/index/star" \
# --rsem_index "${PROJECT_ROOT}/results/genome/rsem" \
# --salmon_index "${PROJECT_ROOT}/results/genome/index/salmon" \
# --transcript_fasta "${PROJECT_ROOT}/results/genome/genome.transcripts.fa" \

printf "Script finished at $(date).\n"
