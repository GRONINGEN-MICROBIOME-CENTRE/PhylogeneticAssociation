#!/bin/bash
#SBATCH --job-name=assign_sgb
#SBATCH --error=phylophlan_assign_sgbs.err
#SBATCH --output=phylophlan_assign_sgbs.out
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10 

module load Anaconda3/2024.02-1
conda activate  /scratch/p300317/metaphlan_sergio/

phylophlan_metagenomic -i Genomes -d SGB.Jan21 -n 1 --nproc 10 -d SGB.Jan19 -o Assignment/ -e fasta
