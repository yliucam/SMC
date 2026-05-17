#!/bin/bash

#SBATCH -A uoa03349           # Project_account (CHANGE TO YOUR ONE)

#SBATCH --mem=300000M #400000M
#SBATCH --time=15:00:00        # Walltime
#SBATCH --ntasks=6            # Number of tasks
#SBATCH --cpus-per-task=20     # Number of OpenMP threads


#SBATCH --mail-user=yliu894@aucklanduni.ac.nz  # Send notifications to this email (CHANGE TO YOUR ONE)
#SBATCH --mail-type=ALL                        # ALL means receive email when code starts and ends

module load JAGS/4.3.1-gimkl-2022a-mt
module load R/4.3.2-foss-2023a  # Call R
srun Rscript test_sub3_repN1.R             # Run R script