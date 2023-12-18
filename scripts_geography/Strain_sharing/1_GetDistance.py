from pathlib import Path
from subprocess import call
import time

N=0
SGB_do =  []
for SGB in Path("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Make_Associations/Association/Results/Geography_cor").glob("*.tsv"):
	SGB_do.append(SGB.stem)

Thresholds = {}
Missing = ""
with open("Data/VallesColomerM_2022_Jan21_thresholds.tsv") as I:
	for line in I:
		if N == 0: N+=1; continue
		l = line.rstrip().split()
		Thresholds[l[0]] = l[1] 
def Extract_distance(SGB, Tree, Send_job=True, Queue="CIBIO"):
	script = "scripts/Extract_norm_distance.R {Tree} {SGB}".format(Tree =str(Tree), SGB=str(SGB))
	
	for_job = "qsub -q {Queue}_cpuQ -o logs/{SGB}_threshold.o -e logs/{SGB}_threshold.e -l mem=10GB -l ncpus=1 -v Tree={Tree},SGB={SGB} scripts/Submit_Distance.sh".format(Tree =str(Tree), SGB=str(SGB), Queue=str(Queue))
	print(for_job)
	print(script)
	if Send_job == True:
		call(for_job, shell=True)

def Complete_missing():
	N = 0
	Q = "CIBIOCM"
	with open("Missing_thresholds.txt") as F:
		for SGB in F:
			SGB = SGB.rstrip()
			Out = "Thresholds/"+SGB+".tsv" 	
			if Path(Out).exists(): continue
			print(SGB)
			Tree = list(Path("../Symlink_phylo").glob("*"+SGB+".TreeShrink.tre"))[0]
			Tree = Tree.resolve()
			for_job = "qsub -q {Queue}_cpuQ -o logs/{SGB}_distance.o -e logs/{SGB}_distance.e -l mem=10GB -l ncpus=1 -v Tree={Tree},SGB={SGB} scripts/Submit_Treshold.sh".format(Tree =str(Tree), SGB=str(SGB), Queue=str(Q))
			call(for_job, shell=True)
			N += 1
			if N == 30: Q = "short"
			if N == 60: Q = "common"
			if N == 90: Q = "CIBIO"
			if N == 120:
				print("Waiting 5 min to submit more jobs")
				time.sleep(300)
				Q = "CIBIOCM"
				N=0

Complete_missing()
exit()
N = 0
Q = "CIBIOCM"
for Tree in Path("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Symlink_phylo/").glob("*.tre"):
	SGB = Tree.stem.split(".")[1]
	if SGB not in SGB_do: continue

	if SGB not in Thresholds: Missing += SGB + "\n"

	if Path("Normalized_distances/"+SGB+".tsv").exists(): pass
	else:  Extract_distance(SGB, Tree, Send_job=True, Queue=Q)


	N += 1
	if N == 30: Q = "short"
	if N == 60: Q = "common"
	if N == 90: Q = "CIBIO"
	if N == 120:
		print("Waiting 5 min to submit more jobs")
		time.sleep(300)
		Q = "CIBIOCM"
		N=0 



with open("Missing_thresholds.txt", "w") as O:
	O.write(Missing)
	
		
