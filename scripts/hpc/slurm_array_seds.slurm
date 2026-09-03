#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 23:00:00
#SBATCH -o logs/job_sed_%a.out
#SBATCH -e logs/job_sed_%a.out
#SBATCH --array=1-32
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --mail-type=END

OFFSET=$1
INDEX=$((SLURM_ARRAY_TASK_ID + OFFSET))

#module load singularity/4.1.0-slurm
singularity exec -B ~/scratch/:/workdir prospectlrdlinux-r.sif \
        Rscript /workdir/compute_fit_seds.R $INDEX '/workdir' 'photoZ_bc03_snorm_LRDMBB_norm_eazy_lrd'
singularity exec -B ~/scratch/:/workdir prospectlinux-r.sif \
        Rscript /workdir/compute_fit_seds.R $INDEX '/workdir' 'photoZ_bc03_snorm_fritz_norm_eazy_lrd'
#singularity exec -B ~/scratch/:/workdir prospectlinux-r.sif \
#        Rscript /workdir/compute_fit_seds.R $INDEX '/workdir' 'photoZ_bc03_snorm_fritz_norm_eazy'
#singularity exec -B ~/scratch/:/workdir prospectlinux-r.sif \
#       Rscript /workdir/compute_fit_seds.R $INDEX '/workdir' 'photoZ_bc03_snorm_norm_eazy'
