from pathlib import Path
import bz2
import pickle

with bz2.open('/shares/CIBIO-Storage/CM/scratch/databases/metaphlansgb_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl', 'r') as rf:
    db = pickle.load(rf)
#tax aqui es la full taxonomy (del k__ al t__) del sgb
dic_size ={}
for tax in db['taxonomy']:
    avg_length = db['taxonomy'][tax][1]
    dic_size[tax.split("|")[-1] ]  = str(avg_length)
with open("Average_genome_length.tsv", "w") as F:
	F.write("SGB\tavg_nt_length\n")
	for tax in dic_size:
		F.write(tax + "\t" + dic_size[tax] + "\n")



