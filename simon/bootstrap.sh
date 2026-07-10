#!/bin/bash

#SBATCH -J dc_boot-job
#SBATCH -c 15
#SBATCH -o dc_boot.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=k.bratkova@lancaster.ac.uk
#SBATCH --mem 50000

srun Rscript /beegfs/client/default/bratkova/STOR_603_Notes/simon/final_iterative_model.R
