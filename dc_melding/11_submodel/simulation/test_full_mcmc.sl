#!/bin/bash

#SBATCH -A uoa03349           # Project_account (CHANGE TO YOUR ONE)

#SBATCH --mem=50000M
#SBATCH --time=0:50:00        # Walltime
#SBATCH --ntasks=100            # Number of tasks
#SBATCH --cpus-per-task=1     # Number of OpenMP threads


#SBATCH --mail-user=yliu894@aucklanduni.ac.nz  # Send notifications to this email (CHANGE TO YOUR ONE)
#SBATCH --mail-type=ALL                        # ALL means receive email when code starts and ends

module load JAGS/4.3.1-gimkl-2022a-mt
module load R/4.3.2-foss-2023a  # Call R
srun Rscript test_full_mcmc.R            # Run R script