#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 2-05:00:00
#SBATCH -o logs/job_lrd_%a.out
#SBATCH -e logs/job_lrd_%a.out
#SBATCH --array=1-55
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --mail-user=jordanacc001@gmail.com
#SBATCH --mail-type=END

OFFSET=$1
INDEX=$((SLURM_ARRAY_TASK_ID + OFFSET))

module load singularity/4.1.0-slurm
singularity exec -B ~/scratch/:/workdir prospectlrdlinux-r.sif \
        Rscript /workdir/prospect_LRD_run.R $INDEX '/workdir' 'photoZ_bc03_snorm_LRDMBB_norm_eazy' 1:55
