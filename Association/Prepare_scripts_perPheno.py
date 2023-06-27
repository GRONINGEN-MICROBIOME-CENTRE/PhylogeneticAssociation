from pathlib import Path
from subprocess import call


def Send_jobs_SGB(SGB):
	#1. Find finished: Everything with a plot
	Done = []
	for File in Path("Results/Plots/"+SGB).glob("*"):
		Done.append(File.stem)
	#2. Get phylogeny
	Phylo = "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Pipeline/Phylogenies/Symlink_phylo/RAxML_bestTree.{SGB}.TreeShrink.tre".format(SGB=SGB)
	#3. For each phenotype that is meant to be run per SGB
	with open( "Results/List_phenos/"+SGB) as F:
		for Pheno in F:
			Pheno = Pheno.rstrip()
			if Pheno == "value": continue
			if Pheno in Done: continue
			#Script text
			Command = """#!/bin/bash
#SBATCH --job-name={SGB}_{Pheno}.job
#SBATCH --output=logs/{SGB}_{Pheno}.o
#SBATCH --error=logs/{SGB}_{Pheno}.e
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=4

ml Anaconda3/2022.05 
source activate /groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Pipeline/Association/R_anpan
ml tbb/2018_U5-GCCcore-7.3.0

Rscript Association_paral.R {SGB} {Tree} {Pheno}
""".format(SGB=SGB, Tree=Phylo, Pheno=Pheno)
			#Name script
			Script_name = "Jobs/{S}_{P}.sh".format(S=SGB, P=Pheno)
			#Create script
			with open(Script_name, "w") as O:
				O.write(Command)
			#4. Send script
			do_command = "sbatch "+Script_name
			print(do_command)
			#call(do_command, shell=True)


for SGB in Path("Results/List_phenos/").glob("*"):
	Send_jobs_SGB(SGB.name)




