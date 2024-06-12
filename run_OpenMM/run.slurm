#!/bin/bash
#SBATCH -c 1					  # Number of Cores per Task
#SBATCH --mem=8192 				  # Requested Memory
#SBATCH -t 240:00:00
#SBATCH -p ials-gpu 					  # Partition
#SBATCH -G 1 					  # Number of GPUs
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --job-name=wt-300K

python -u run.nvt.py *.pdb *.psf 0		
