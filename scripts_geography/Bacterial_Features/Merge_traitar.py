from glob import glob, iglob
import pandas as pd
import numpy as np
SGBs_to_run = []
with open("../Association/Tree_sizes2.txt") as F:
        for line in F:
                l = line.split()
                if l[0] == "Tree": continue
                if int(l[1]) < 300: continue
                SGBs_to_run.append(l[0].split(".")[1])

SGBs_done = []
with open("../Phenotypes/traitar_output_646SGBs_tab.txt") as F:
        for line in F:
                l = line.split()[0]
                if l == "SGB": continue
                SGBs_done.append( l)

selected_sgbs = list(set(SGBs_to_run) - set(SGBs_done))



merged_sgbs = list()
with open('/shares/CIBIO-Storage/CM/scratch/databases/chocophlansgb/releases/Jan21/Jan21_merged_sgbs.csv', 'r') as read_file:
    read_file.readline()
    for line in read_file:
        line = line.strip().split(';')
        merged_sgbs.append([line[0][3:]] + line[1].split(','))

dfs = []
for sgb in selected_sgbs:
    sgb_id = sgb.replace('t__SGB','').replace('_group','')
    found = False
    sgbs = [sgb_id]
    for merge in merged_sgbs:
        if sgb_id in merge:
            sgbs = merge
            
    for folder in iglob('Traitar/*'):
        name = folder.split('/')[-1]
        if name in sgbs:
            found = True
            file = folder +'/traitar/phenotype_prediction/predictions_flat_majority-votes_combined.txt'
            break
    if not found:       
        for folder in iglob('/shares/CIBIO-Storage/CM/scratch/projects/nkarcher_FMT_meta/results/traitar/output/*'):
            name = folder.split('/')[-1]
            if name in sgbs:
                found = True
                file = folder +'/traitar/phenotype_prediction/predictions_flat_majority-votes_combined.txt'
                break
    if not found:       
        for folder in iglob('/shares/CIBIO-Storage/CM/scratch/users/aitor.blancomiguez/tests/traitar/results/*'):
            name = folder.split('/')[-1]
            if name in sgbs and os.path.exists(folder +'/traitar/phenotype_prediction/predictions_flat_majority-votes_combined.txt'):
                found = True
                file = folder +'/traitar/phenotype_prediction/predictions_flat_majority-votes_combined.txt'
                break
                
    if found:
        df = pd.read_csv(file, sep='\t', names=['sample', 'phenotype', 'value', 'model'], skiprows=1)
        df['sgb'] = sgb
        del df['sample']
        dfs.append(df)
    else: print(sgb_id)

df = pd.concat(dfs)
df = df.groupby(['sgb', 'phenotype']).sum()['value'] / 2
df_traitar = df.reset_index().pivot(index='sgb', columns='phenotype', values='value').fillna(0).replace(0.5, np.nan)


#Merge with Mireia's
#M = pd.read_csv("../Phenotypes/traitar_output_646SGBs_tab.txt", sep="\t")
#df_traitar = pd.concat([df_traitar, M], ignore_index=False)


df_traitar.to_csv("traitar_output2.txt", sep="\t",na_rep="NA" )
