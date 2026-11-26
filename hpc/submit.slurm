#!/bin/bash
# SLURM array job for the full simulation. Adjust partition, account, time and
# memory to your cluster. Run `Rscript hpc/run_task.R calibrate` first; it prints
# the --array range. Blocks differ a lot in cost, so you may prefer to submit
# Block C tasks separately with more memory and time (see results/tasks.csv).
#SBATCH --job-name=survsim
#SBATCH --array=1-1000%200          # replace 1000 with the number of tasks; %200 = max concurrent
#SBATCH --cpus-per-task=4           # 4 MCMC chains in parallel
#SBATCH --mem=16G                   # Block C (n_test = 50,000, p = 2,003) needs about 16 GB
#SBATCH --time=48:00:00
#SBATCH --output=logs/survsim_%A_%a.out

module load R                       # site-specific
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$SLURM_SUBMIT_DIR"
Rscript hpc/run_task.R task "$SLURM_ARRAY_TASK_ID"
