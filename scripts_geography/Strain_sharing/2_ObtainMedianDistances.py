from pathlib import Path
from subprocess import call
import time

def Send_job(SGB, FILE, Send_job=True, Queue="CIBIO"):
        script = "scripts/Compare_distances_geography.R {SGB}".format(SGB=str(FILE))
        for_job = "qsub -q {Queue}_cpuQ -o logs/{SGB}_median.o -e logs/{SGB}_median.e -l mem=10GB -l ncpus=1 -v SGB={File} scripts/Submit_Median.sh".format(SGB=str(SGB), Queue=str(Queue), File=str(FILE))
        print(for_job)
        print(script)
        if Send_job == True:
                call(for_job, shell=True)

Q = "CIBIOCM"
N=0
for Distance_file in Path("Normalized_distances").glob("*.tsv"):
	SGB = Distance_file.stem
	Output = "Distance_countryNcontinent/"+SGB+"_continent.tsv"
	if Path(Output).exists(): continue
	Send_job(SGB, Distance_file, Send_job=True, Queue=Q)
	N += 1
	if N == 30: Q = "short"
	if N == 60: Q = "common"
	if N == 90: Q = "CIBIO"
	if N == 120:
		print("Waiting 2 min to submit more jobs")
		time.sleep(120)
		Q = "CIBIOCM"
		N=0 

