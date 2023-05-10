import pickle
import re
import bz2
from pathlib import Path
import numpy as np

Hash_size = "/shares/CIBIO-Storage/CM/cmstore/tools/anaconda3/envs/metaphlan-4/lib/python3.9/site-packages/MetaPhlAn4beta1-4b1-py3.9.egg/metaphlan4beta1/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl"
with bz2.open(Hash_size, 'rb') as f:
        data = pickle.load(f, fix_imports=True)
Hash_size = data["taxonomy"]

Reads_per_sample = {}
Counts = "Metadata_fecal_selection.tsv"
counter_l = 0
with open(Counts) as f:
    for line in f:
        if counter_l == 0: counter_l += 1 ; continue
        l = line.rstrip().split()
        try: Reads_per_sample[l[1]] = [l[1], float(l[16]), float(l[19]) ]
        except: continue
        #Sample: cohort, read number, median read length    

def Read(Already_made = False ):
    Samples_made = []
    if Already_made == True:
        with open("PublicDepth.tsv") as F:
                for line in F:
                    l = line.rstrip().split()
                    Samples_made.append(l[0])
        Samples_made = set(Samples_made)                    

    Missing_sample = []
    Dic_depth = {}
    with open("File_location.txt") as F:
        for line in F:
            l = line.rstrip().split()
            if Already_made == True:
                if Sample in Samples_made: continue
            Sample = l[0]
            if Sample not in Reads_per_sample: Missing_sample.append(Sample)
            Dic_depth[Sample] = {}
            Cohort = l[1]
            MP4 = l[2]
            SF4 = l[3]
            if not Path(MP4).exists(): continue
            with open(MP4) as F2:
                for line2 in F2:
                    if line2[0] == "#": continue
                    l2 = line2.rstrip().split()
                    Taxa = l2[0]
                    if not "|t__" in Taxa: continue
                    Reads_clade = float(l2[2])
                    
                    Genome_size =  float(Hash_size[Taxa][1])
                    Info = Reads_per_sample[Sample]
                    Number_reads = float(Info[1])
                    Read_length = float(Info[2])

                    Depth_coverage = Read_length*(Number_reads*Reads_clade) / Genome_size
                    Dic_depth[Sample][Taxa] = str(Depth_coverage)


def Save_file(Dic_depth, Output_name):
    print("Saving")
    to_write = ""
    for key in Dic_depth:
        for sgb in Dic_depth[key]:
            coverage_sgb = Dic_depth[key][sgb]
            if coverage_sgb == 0 : continue
            to_write += "{Sample}\t{Bug}\t{Abundance}\n".format(Sample=key, Bug=sgb, Abundance=str(coverage_sgb))
        with open(Output_name, "w") as O:
            O.write(to_write)
def Attach_file(Dic_depth, Output_name):
    print("Saving")
    to_write = ""
    for key in Dic_depth:
        for sgb in Dic_depth[key]:
            coverage_sgb = Dic_depth[key][sgb]
            if coverage_sgb == 0 : continue
            to_write += "{Sample}\t{Bug}\t{Abundance}\n".format(Sample=key, Bug=sgb, Abundance=str(coverage_sgb))
            with open(Output_name, "a") as O:
                O.write(to_write)

Save_file(Dic_depth, "PublicDepth.tsv")
print("Samples missing from metadata, possibly due to unknown number of reads or read length: " + " ".join(Missing_sample) )
