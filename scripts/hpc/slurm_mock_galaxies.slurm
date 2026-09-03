#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -t 03:00:00
#SBATCH -o logs/job_mock.out
#SBATCH -e logs/job_mock.out
#SBATCH --ntasks=32
#SBATCH --mem=20G

module load singularity/4.1.0-slurm

singularity exec -B ~/scratch/:/workdir prospectlinux-r.sif \
        Rscript /workdir/mock_galaxies_uniform.R 

