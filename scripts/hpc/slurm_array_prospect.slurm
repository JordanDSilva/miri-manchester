#!/bin/bash --login
#SBATCH -p work
#SBATCH -t 09:00:00
#SBATCH -o logs/job_%a.out
#SBATCH -e logs/job_%a.out
#SBATCH --array=1-800
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --mail-user=jordanacc001@gmail.com
#SBATCH --mail-type=END

OFFSET=$1
INDEX=$((SLURM_ARRAY_TASK_ID + OFFSET))

module load singularity/4.1.0-slurm
singularity exec -B /scratch/pawsey0335/jdsilvauwa/:/workdir prospectlinux-r.sif \
        Rscript /workdir/prospect_run.R $INDEX '/workdir' 'photoZ_bc03_snorm' 
