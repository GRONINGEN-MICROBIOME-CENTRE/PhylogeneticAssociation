import time
from pathlib import Path
from subprocess import call

script_file = "/shares/CIBIO-Storage/CM/scratch/users/aitor.blancomiguez/tools/strain_transmission/strain_transmission_stats.py"
script2 = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Strain_sharedness/scripts/Get_proportion_sharedness.py"
Thresholds_file = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Strain_sharedness/Data/Jan21_thresholds.tsv"
Tree_dir = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Symlink_phylo"
Metadata = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Strain_sharedness/Data/Input_intrapopulation.tsv"
Dist_dir = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Strain_sharedness/Normalized_distances"

Env = "/shares/CIBIO-Storage/CM/scratch/tools/anaconda3/envs/metaphlan_sergio"


Q = "CIBIO"

Info = """
Start executionusage: strain_transmission.py [-h] [-t TREE] [-d DIST] [-m METADATA] [-o OUTPUT_DIR] [--save_dist] [--csv] [--vertical] [--horizontal] [--restrictive] [--permissive] [--threshold THRESHOLD] [--projects PROJECTS [PROJECTS ...]]

options:
  -h, --help            show this help message and exit
  -t TREE, --tree TREE  The input tree file
  -d DIST, --dist DIST  The input PhyPhlAn pairwise distances file
  -m METADATA, --metadata METADATA
                        The input metadata
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        The output directory
  --save_dist           [Optional] Save the PhyPhlAn pairwise distances file
  --csv                 [Optional] Export the results as CSV
  --vertical            [Optional] Detect only vertical transmission events
  --horizontal          [Optional] Detect only horizontal transmission events
  --restrictive         [Optional] Use only one timepoint of one individual per family/house in the transmission threshold selection
  --permissive          [Optional] Use all the data in the transmission threshold selection
  --threshold THRESHOLD
                        [Optional] A custom distribution threshold value. Default: 0.01
  --projects PROJECTS [PROJECTS ...]
                        [Optional] The specific projects in which detect the transmission

"""



#Get path tree per SGB
Tree_path = {}
for Tree in Path(Tree_dir).glob("*.tre"):
    SGB = Tree.stem.split(".")[1]
    Tree_path[SGB] = Tree.resolve()
#Get threshold per SGB
Thresholds = {}
with open(Thresholds_file) as F:
    for line in F:
        l = line.strip().split()
        Thresholds[l[0]] = l[1]



#Iterate distances
N = 0
Commands = ""
for Distance_file in Path(Dist_dir).glob("*.tsv"):
    SGB = Distance_file.stem
    Tree_sgb = Tree_path[SGB]
    if SGB not in Thresholds: Threshold = 0.03 #continue
    else: Threshold = Thresholds[SGB]
    
    Output_dir = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Strain_sharedness/Transmission/"+SGB
    if Path(Output_dir + "/Tranmission_proportions.tsv").exists(): pass #continue
    if not Path(Output_dir).exists(): Path(Output_dir).mkdir(parents=True, exist_ok=True)
    Command = "python {strain_script} -d {dist} -m {meta} -o {out} --threshold {thresh}".format(strain_script=script_file, Tree=Tree_sgb,  dist=Distance_file, meta=Metadata, out=Output_dir, thresh=Threshold)
    Command2 = "python  {script} {SGB}".format(script=script2, SGB=SGB)
    PBS_command = """#!/bin/bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/shares/CIBIO-Storage/CM/scratch/tools/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
        eval "$__conda_setup"
else
        if [ -f "/shares/CIBIO-Storage/CM/scratch/tools/anaconda3/etc/profile.d/conda.sh" ]; then 
                . "/shares/CIBIO-Storage/CM/scratch/tools/anaconda3/etc/profile.d/conda.sh"
        else
                export PATH="/shares/CIBIO-Storage/CM/scratch/tools/anaconda3/bin:$PATH"
        fi
fi
unset __conda_setup


# <<< conda initialize <<<
conda deactivate
conda activate {Env}

cd /shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Strain_sharedness
""".format(Env = Env) + "\n" + Command + "\n" + Command2
    
    Script = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Strain_sharedness/scripts/tmp/{SGB}.sh".format(SGB=SGB)
    with open(Script, "w") as S:
        S.write(PBS_command)
    Send_command = "qsub -q {Queue}_cpuQ -o logs/{SGB}_share.o -e logs/{SGB}_share.e -l mem=10GB -l ncpus=1 {script}".format(SGB=str(SGB), Queue=str(Q), script=Script )
    Commands += Send_command +"\n"
    print(Commands)
    call(Send_command, shell=True)

    continue    
    N += 1
    if N == 30: Q = "short"
    if N == 60: Q = "common"
    if N == 90: Q = "CIBIO"
    if N == 120:
        print("Waiting 5 min to submit more jobs")
        time.sleep(300)
        Q = "CIBIOCM"
        N=0 


