import pickle
import re
import bz2

Hash_size = "/groups/umcg-tifn/tmp01/tools/conda_biobakery4/lib/python3.1/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl"
with bz2.open(Hash_size, 'rb') as f:
	data = pickle.load(f, fix_imports=True)
Hash_size = data["taxonomy"]

Reads_per_sample = {}
Counts = "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Pipeline/Files/Read_count/Merged_counts.tsv"
with open(Counts) as f:
	for line in f:
		l = line.rstrip().split()
		Reads_per_sample[l[0]] = l[1]


#All cohorts
Cohorts= "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Pipeline/Files/Merged_Groningen.tsv"
#LLnext
Next = "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Pipeline/Files/Next/LLNEXT_metaphlan_4_CLEAN_03_03_2023.txt"

def Calculate_depth(Hash_size, Reads_per_sample, Tax_abu):
	Dic_depth = {}
	Position = 0
	Taxa_index = {}
	Missing = []
	with open(Tax_abu) as f:
		for line in f:
			l = line.rstrip().split()
			if Position == 0:
				for i in range(len(l)):
					Taxa= l[i]
					if ".t__" in Taxa: Taxa = re.sub(r'(\w+)\.(?=\w__)', r'\1|' , Taxa)
					Taxa_index[i] = Taxa
			else:
				for i in range(len(l)):
					Taxa = Taxa_index[i].strip('"')
					if Taxa == "ID": 
						Sample = l[i].strip('"').replace("_metaphlan", "").replace("_FU", "_F").replace("_APK", "")
						if Sample[0:2] == "HV" and "FSK" in Sample: Sample = Sample.split("_")[0]
						if Sample not in Reads_per_sample: Missing.append(Sample) ; break
						print(Sample)
						Number_reads = int(Reads_per_sample[Sample])
						Dic_depth[Sample] = {}
					elif Taxa == "UNCLASSIFIED": 
						Percentage_unknown = float(l[i])
						#Number_reads_map = Number_reads*Percentage_unknown/100
					else:
						#if Sample in Missing: continue
						if not "t__" in Taxa or  "k__Eukaryota" in Taxa: continue
						Genome_size =  float(Hash_size[Taxa][1])
						Reads_clade =  float(l[i])
						Depth_coverage = 150*(Number_reads*Reads_clade) / Genome_size
						#if Reads_clade>0 :  print(Genome_size, l[i],Taxa)
						Dic_depth[Sample][Taxa] = Depth_coverage
			Position += 1
	print("Missing samples (in MP4 but not in number of sample files) : "+ " ".join(Missing))
	return(Dic_depth, Missing)

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
	


#NEXT
Next_quan,Missing_next = Calculate_depth(Hash_size, Reads_per_sample, Next)
Save_file(Next_quan, "Baby_DepthCoverage.tsv")

#Rest
Quan,Missing = Calculate_depth(Hash_size, Reads_per_sample, Cohorts)
Save_file(Quan, "DepthCoverage.tsv")

					
					
			











