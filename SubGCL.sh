#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --ntasks=16
#SBATCH --mem=160g
#SBATCH --tmp=50g
#SBATCH --mail-type=BEGIN,END,FAIL  
#SBATCH --mail-user=bartz209@umn.edu

module load conda

# Initialize Conda for this shell
conda activate PaolaR-env

cell_type=$1
sample_id=$2
Rscript GCLRunner.R "$cell_type" "$sample_id"
