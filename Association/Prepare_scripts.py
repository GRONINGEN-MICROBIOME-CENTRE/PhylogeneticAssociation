from pathlib import Path
from subprocess import call


#1. Find finished
Done = []
for File in Path("Results/Reports/").glob("*"):
	Done.append(File.stem)


#2. Find trees
# /groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Pipeline/Phylogenies/Symlink_phylo/
Phylo_dir = "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Pipeline/Phylogenies"
for Tree in Path(Phylo_dir).rglob("*.tre"):
	SGB = "t__" + Tree.stem.split("t__")[1].split(".")[0]
	if SGB in Done: continue
	Command = """#!/bin/bash
#SBATCH --job-name={SGB}.job
#SBATCH --output=logs/{SGB}.o
#SBATCH --error=logs/{SGB}.e
#SBATCH --time=16:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=4

ml Anaconda3/2022.05 
source activate /groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Pipeline/Association/R_anpan
ml tbb/2018_U5-GCCcore-7.3.0

Rscript Association.R {SGB} {Tree}

""".format(SGB=SGB, Tree=str(Tree))
	Script_name = "Jobs/{S}.sh".format(S=SGB)
	with open(Script_name, "w") as O:
		O.write(Command)
	do_command = "sbatch "+Script_name
	
	Command2 = """#ml Anaconda3/2022.05 
#source activate /groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Pipeline/Association/R_anpan
#ml tbb/2018_U5-GCCcore-7.3.0

Rscript 0_List_phenotypes_per_SGB.R  {SGB} {Tree}""".format(SGB=SGB, Tree=str(Tree))
	#Run command to list phenotypes
	print(Command2)
	call(Command2, shell=True)	
	#Run command for script
	#call(do_command, shell=True)




