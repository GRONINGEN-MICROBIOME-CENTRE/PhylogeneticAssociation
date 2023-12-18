from pathlib import Path
import pandas as pd
from io import StringIO
from subprocess import call

def Get_info(SGB):
	Info = {}
	with open("logs/{SGB}.o".format(SGB=SGB)) as O:
		for line in O:
			#[1] "Age for association"
			if "for association" in line: pheno = line.split()[1].strip('"') ; Info[pheno] = {"ok":"0", "bad":"0", "very_bad":"0" }
			elif "(ok)" in line:  N = line.split()[3]; Info[pheno]["ok"] = N
			elif "(bad)" in line:  N = line.split()[3]; Info[pheno]["bad"] = N
			elif "(very bad)" in line:  N = line.split()[4]; Info[pheno]["very_bad"] = N
	return(Info)

def Check_if_exist(File, N_param):
	if File.exists():
		n=0
		with open(File) as I:
			for line in I:
				if n==0: n+=1 ; continue
				P = line.rstrip().split()
	else: 
		P = ["NA" for _ in range(N_param)]
	return(P)

def Merge(Association_folder, Output_file, Print=False):
	Phenos_SGB = []
	if Path(Output_file).exists():
		N = 0
		Associations = ""
		with open(Output_file) as F:
			for line in F:
				Associations += line.replace(",", "\t")
				l = line.split(",")
				Phenos_SGB.append( l[4] + "__" + l[5])
	else:
		Associations = "elpd_diff\tse_diff\tAbs_LowerBound\tT_stat\tPhenotype\tSGB\tCovariates\tk_good\tk_ok\tk_bad\tk_verybad"
		Phylogenetic_effects = "\t" + "\t".join(["phylo_mean","phylo_median","phylo_sd","phylo_mad","phylo_q5","phylo_q95","phylo_rhat","phylo_ess_bulk","phylo_ess_tail"]) + "\n"
		Associations = Associations + Phylogenetic_effects
	
	for File in Path("Results/"+Association_folder+"/Reports/").glob("*.tsv"):
		with open(File) as F:
			if Print == True: print(File.name)
			N = 0
			for line in F:
				if N == 0: N+=1; continue
				l = line.rstrip().split("\t")
				SGB = l[-1]
				Phenotype = l[-3]
				if Phenotype + "__" + SGB in Phenos_SGB: continue
				Pareto_scores = Path("Results/{subfolder}/Models/{SGB}/{Phenotype}/pareto.tsv".format(subfolder=Association_folder, SGB=SGB, Phenotype=Phenotype))
				Pareto_scores = Check_if_exist(Pareto_scores, 4)
				Phylo_Estimates = Path("Results/{subfolder}/Models/{SGB}/{Phenotype}/phylo_estimates.tsv".format(subfolder=Association_folder, SGB=SGB, Phenotype=Phenotype))
				Phylo_Estimates = Check_if_exist(Phylo_Estimates, 10)		

				Associations_t = [l[0], l[1], str(abs(float(l[0])) - 2*abs(float(l[1]))), str(abs(float(l[0])/float(l[1]))), Phenotype, SGB, l[-2]]
				Associations_t.extend(Pareto_scores)
				Associations_t.extend(Phylo_Estimates[1:])
				Associations += "\t".join(Associations_t) + "\n"
				
			

	data_file = StringIO(Associations)
	df = pd.read_csv(data_file, delimiter='\t')
	sorted_df = df.sort_values('T_stat', ascending=False)
	sorted_df.to_csv(Output_file, index=False)

	print("Filtering results")
	Get_significant = "Rscript Useful_scripts/Explore_Assocations.R "+Output_file.replace(".csv","")
	call(Get_significant, shell=True)

Merging_dic = { "Age,Continent" : 'All_results.csv', "Age": 'Geography_results.csv', 'Comparison_SequencingProtocol': "SequencingProtocol_results.csv", "Disease_specific/Age,Continent" : "Results_disease.csv", "Disease_specific/Age,Continent,BMI" : "Disease_specific/Results_hypertension_BMI.csv", "Disease_specific/Age,BMI_Asia":"Results_hypertension_BMI_Asia.csv", "Disease_specific/Age,BMI_Europe":"Results_hypertension_BMI_Europe.csv", "Disease_specific/Age_Asia": "Results_disease_Asia.csv", "Disease_specific/Age_Europe": "Results_disease_Europe.csv", "Disease_specific/Age_North_America": "Results_disease_NorthAmerica.csv", "Disease_per_Continent/Age,Country,BMI_Asia": "Results_PerContinent_hypetension_BMI_Asia.csv", "Disease_per_Continent/Age,Country,BMI_Europe":"Results_PerContinent_hypetension_BMI_Europe.csv", "Disease_per_Continent/Age,Country_Asia":"Results_PerContinent_Asia.csv","Disease_per_Continent/Age,Country_Europe":"Results_PerContinent_Europe.csv", "Disease_per_Continent/Age,Country_North_America":"Results_PerContinent_NorthAmerica.csv" }
#Merging_dic = { "Disease_specific/Age_Asia": "Results_disease_Asia.csv", "Disease_specific/Age_Europe": "Results_disease_Europe.csv", "Disease_specific/Age_North_America": "Results_disease_NorthAmerica.csv" }

#Merging_dic = { "Disease_specific/Age,Continent" : "Results_disease.csv", "Disease_specific/Age,Continent,BMI" : "Results_hypertension_BMI.csv", "Disease_specific/Age,BMI_Asia":"Results_hypertension_BMI_Asia.csv", "Disease_specific/Age,BMI_Europe":"Results_hypertension_BMI_Europe.csv", "Disease_specific/Age_Asia": "Results_disease_Asia.csv", "Disease_specific/Age_Europe": "Results_disease_Europe.csv", "Disease_specific/Age_North_America": "Results_disease_NorthAmerica.csv", "Disease_per_Continent/Age,Country,BMI_Asia": "Results_PerContinent_hypetension_BMI_Asia.csv", "Disease_per_Continent/Age,Country,BMI_Europe":"Results_PerContinent_hypetension_BMI_Europe.csv", "Disease_per_Continent/Age,Country_Asia":"Results_PerContinent_Asia.csv","Disease_per_Continent/Age,Country_Europe":"Results_PerContinent_Europe.csv", "Disease_per_Continent/Age,Country_North_America":"Results_PerContinent_NorthAmerica.csv" }

#Merging_dic = {"Continent_stratified/Age,Country_Asia":"Results_PerContinent_Asia.csv", "Continent_stratified/Age,Country_Europe":"Results_PerContinent_Europe.csv", "Continent_stratified/Age,Country_North_America":"Results_PerContinent_NorthAmerica.csv" }

#Merging_dic = {"Continent_stratified/Age,BMI_Europe/": "Results_PerContinent_BMI_Europe.csv", "Continent_stratified/Age,BMI_Asia/": "Results_PerContinent_BMI_Asia.csv" }
#Merging_dic = {"Disease_per_Continent/Age,Country,BMI_Europe" :  "Results_Disease_BMI_Europe.csv" }

#Merging_dic = {"Continent_stratified/Age,Country_Europe/": "Results_PerContinent_Europe.csv", "Continent_stratified/Age,Country_Asia/": "Results_PerContinent_Asia.csv", "Continent_stratified/Age,Country_Africa/" : "Results_PerContinent_Africa.csv", "Continent_stratified/Age,Country_North_America/" : "Results_PerContinent_NorthAmerica.csv",  "Continent_stratified/Age,Country_South_America/" : "Results_PerContinent_SouthAmerica.csv", "Continent_stratified/Age,Country_Oceania/" : "Results_PerContinent_Oceania.csv"   }
#Merging_dic = {"Abundance_association/": "Results_Abundance.csv"}
#Merging_dic = {"Cohort_specific/Age_XuQ_2021/":"Results_XuQ_Age.csv"}
#Merging_dic = {"Continent_stratified/Country_Asia/": "Results_PerContienent_Asia_ElderInfant.csv", "Continent_stratified/Country_Europe/": "Results_PerContienent_Europe_ElderInfant.csv", "Continent_stratified/Country_Africa/": "Results_PerContienent_Africa_ElderInfant.csv", "Continent_stratified/Country_North_America/": "Results_PerContienent_NorthAmerica_ElderInfant.csv" }
Merging_dic = {"Continent_stratified/Age,Country_Europe/UC_CD/": "Results_UC_CD.csv", "Cohort_specific/Age_DAG3,IBD,LLD/": "Results_DutchIBD.csv" }

#Merging_dic = {"Cohort_specific/Age_DAG3,IBD,LLD/": "Results_DutchIBD.csv"}

Do_individual = False

#Do individual one
if Do_individual == True:
	#Association_folder = "Disease_specific/Age,Continent"
	Association_folder = "Age,Continent"
	#Association_folder = "Age"
	#Association_folder = "Comparison_SequencingProtocol"
	Output_file = Merging_dic[Association_folder]
	print("Merging dir {fol} to {o}".format(fol=Association_folder, o=Output_file))
	Merge(Association_folder, Output_file, Print=False)
else:
	for Association_folder in Merging_dic:
		Output_file = Merging_dic[Association_folder]
		print("Merging dir {fol} to {o}".format(fol=Association_folder, o=Output_file))
		Merge(Association_folder, Output_file, Print=False)




