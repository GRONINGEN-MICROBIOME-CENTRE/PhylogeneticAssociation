from pathlib import Path
import pandas as pd
from io import StringIO


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

#Find all reports
Associations = "elpd_diff\tse_diff\tAbs_LowerBound\tT_stat\tPhenotype\tSGB\tk_ok\tk_bad\tk_verybad\n"
Prev_SGB = ""
for File in Path("Results/Reports/").glob("*.tsv"):
	with open(File) as F:
		print(File.name)
		N = 0
		for line in F:
			if N == 0: N+=1; continue
			l = line.rstrip().split()
			SGB = l[-1]
			if SGB != Prev_SGB: Info = Get_info(SGB) ; Prev_SGB = SGB
			Info_pheno = Info[l[-2]]
			Associations += "\t".join([l[0], l[1], str(abs(float(l[0])) - 2*abs(float(l[1]))), str(abs(float(l[0])/float(l[1]))), l[-2], SGB, Info_pheno["ok"], Info_pheno["bad"], Info_pheno["very_bad"] ] ) + "\n"
			

data_file = StringIO(Associations)
df = pd.read_csv(data_file, delimiter='\t')
sorted_df = df.sort_values('T_stat', ascending=False)

sorted_df.to_csv('All_results.csv', index=False)


