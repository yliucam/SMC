#!/bin/bash

#SBATCH -A uoa03349           # Project_account (CHANGE TO YOUR ONE)

#SBATCH --time=2:00:00        # Walltime
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=1     # Number of OpenMP threads
#SBATCH --mem-per-cpu=20000   # Memory per CPU (in MB)

#SBATCH --mail-user=yliu894@aucklanduni.ac.nz  # Send notifications to this email (CHANGE TO YOUR ONE)
#SBATCH --mail-type=ALL                        # ALL means receive email when code starts and ends

module load JAGS/4.3.1-gimkl-2022a-mt
module load R/4.3.2-foss-2023a  # Call R
srun Rscript cap_recap_run.R              # Run R script