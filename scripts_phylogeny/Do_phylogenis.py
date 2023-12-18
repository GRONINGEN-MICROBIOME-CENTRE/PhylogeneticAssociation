from pathlib import Path
from subprocess import call


Only_do_this = []
with open("Send_agin_nomutations.txt") as F:
	N = 0
	for line in F:
		if N == 0: N+=1 ; continue
		Only_do_this.append(line.rstrip().split()[0])

#1. Choose files to use
Prevalence_info = "/shares/CIBIO-Storage/CM/scratch/users/sergio.andreusanchez/Number_non0.tsv"
Samples_do = []
Threshold =  200
Ordered_samples = {}

queue="CIBIO_cpuQ"
#queue="CIBIOCM_cpuQ"
cpu="10"
memory="25"

Phylo_done = []
for Phylo in Path("Phylogenies").glob("*/RAxML_bestTree.*"):
	Phylo_done.append(str(Phylo.parent.name))
Phylo_running = []
#with open("IDs_running_jobs.txt") as F:
#	for line in F:
#		Phylo_running.append(line.rstrip())

with open(Prevalence_info) as P:
    for line in P:
        l = line.rstrip().split()
        if l[0] == "SGB": continue
        SGB = l[0].split("|")[-1]
        if SGB not in Only_do_this: continue
        if not "t__" in SGB: continue
        N_samples = int(l[1])
        #Another option: change requested cpu and memory depending on number of samples above 0
        if N_samples >= Threshold: 
            Samples_do.append(SGB)
            Ordered_samples[SGB] = N_samples
Ordered_samples = dict(sorted(Ordered_samples.items(), key=lambda x: x[1], reverse=True))
print(len(Phylo_done))

print(len(Samples_do))
Batch_n = 0 #update to Batch for next round
Batch = 20000
N = 0 
R = ""
#2. Prepare job and write script
for SGB_dir in Ordered_samples:
    N += 1
    #if N <= Batch_n: continue
    #if N == Batch: break
    SGB_dir = Path("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Samples_tree/SGB") / SGB_dir
#for SGB_dir in Path("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Samples_tree/SGB").glob("*"):
    SGB_name = SGB_dir.name
    if SGB_name in Phylo_done or SGB_name in Phylo_running: N -= 1; print(SGB_name); continue
    #if SGB_name not in Samples_do: continue
    Primary = SGB_dir / "2_x/Primary"
    Secondary = SGB_dir / "2_x/Secondary"
    DIR = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Phylogenies/" + SGB_name
    if not Path(DIR).exists(): call("mkdir "+DIR, shell=True)
    Script = """#!/bin/bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/shares/CIBIO-Storage/CM/scratch/tools/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/shares/CIBIO-Storage/CM/sratch/tools/anaconda3/etc/profile.d/conda.sh" ]; then 
        . "/shares/CIBIO-Storage/CM/scratch/tools/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/shares/CIBIO-Storage/CM/scratch/tools/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup


# <<< conda initialize <<<
conda deactivate
conda activate /shares/CIBIO-Storage/CM/scratch/tools/anaconda3/envs/metaphlan_sergio
strainphlan -o {OUTPUT} -n {CPU}  --sample_with_n_markers 50 --secondary_sample_with_n_markers 50  --sample_with_n_markers_after_filt 33 --marker_in_n_samples 50 --samples {PRIMARY} --secondary_samples {SECONDARY} -c {SGB} -d /shares/CIBIO-Storage/CM/scratch/databases/metaphlansgb_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl --treeshrink --debug 
""".format(OUTPUT=DIR, PRIMARY=Primary, SECONDARY = Secondary, SGB= SGB_name, CPU=cpu)
    Script_name = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Scripts_phylo/script_" + SGB_name + ".sh"
    #if Path(Script_name).exists(): N-=1 ; continue
    with open(Script_name, "w") as Output:
        Output.write(Script)
    #3. Submit script
    #call("qsub  -q {Q} -o {logs_dir}/{sample}.o -e {logs_dir}/{sample}.e  -l mem={M}GB -l ncpus={CPU} {FILE}".format(Q=queue, CPU=cpu, M=memory, FILE= Script_name, sample=SGB_name, logs_dir="/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Scripts_phylo/logs", shell=True ))
    R += "qsub  -q {Q} -o {logs_dir}/{sample}.o -e {logs_dir}/{sample}.e  -l mem={M}GB -l ncpus={CPU} {FILE}\n".format(Q=queue, CPU=cpu, M=memory, FILE= Script_name, sample=SGB_name, logs_dir="/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Scripts_phylo/logs") 

with open("Run_script_trees.sh", "w") as F:
	F.write(R)

