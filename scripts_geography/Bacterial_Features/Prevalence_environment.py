from pathlib import Path


Locations = """Airways
Freshwater
Lab_Mice
NHP
Ocean
Oral_Ancient
Oral_NW
Oral_W
Rummen
Skin
Soil
Stool_Ancient
Stool_NW
Stool_W
Vagina
Wild_Mice
""".split()

Dic_SGB = {}
#Rummen  SGB37550__GGB25580_SGB37550     288     313
with open("/shares/CIBIO-Storage/CM/scratch/projects/ablanco_metaphlansgb/results/new_diversity/prevalences.tsv") as F:
	for line in F:
		l = line.rstrip().split()
		Location = l[0]
		SGB = "t__" + l[1].split("_")[0]
		if SGB not in Dic_SGB:
			Dic_SGB[SGB] = {}
			for L in Locations: Dic_SGB[SGB][L] = [0, 0, 0] 
		prevalence = float(l[2])
		Total = float(l[3])
		Prevalence_bin = 1
		rel_prevalence = prevalence/Total
		Dic_SGB[SGB][Location] = [prevalence, rel_prevalence, Prevalence_bin]

with open("prevalence_SGB_environments.tsv", "w") as F:
	F.write("SGB\tenvironment\tprevalence\trel_prevalence\tprevalence_bin\n")
	for SGB in Dic_SGB:
		# environment prevalence rel_prevalence prevalence_bin
		for environment in Dic_SGB[SGB]:
			info = Dic_SGB[SGB][environment]
			To_write = "\t".join([SGB, environment, str(info[0]), str(info[1]), str(info[2])]) + "\n"
			F.write(To_write)



