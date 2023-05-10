from pathlib import Path


Data_info = "Metadata_fecal_selection.tsv"
Dic_sample_project = {}
N = 0
with open(Data_info) as F:
    for line in F:
        if N == 0: N+=1 ;  continue
        l = line.rstrip().split()
        if l[0] not in Dic_sample_project: Dic_sample_project[l[0]] = []
        Dic_sample_project[l[0]].append(l[1])
    

def Search():
    return()



Rename = { "ThomasAM_2018a" : "ThomasAM_2019_a", "TettAJ_2019_c":"CM_ethiopia", "DavidLA_2015":"LawrenceA_2015", "TettAJ_2019_b":"CM_ghana", "TettAJ_2019_a":"CM_tanzania", "PasolliE_2019":"CM_madagascar", 'FerrettiP_2018':'CM_caritro', 'ThomasAM_2018b': 'ThomasAM_2018_b', 'Bengtsson-PalmeJ_2015': 'BengtssonPalmeJ_2015', 'GhensiP_2019': 'CM_periimplantitis'}
Additional_rename = {"BaruchEN_2021": "BaruchE_2021", "Valles-ColomerM_2023_a": "VallesColomerM_2023_a", "ThomasAM_2018b":"ThomasAM_2019_b", "LeeKA_2022":"CM_SEERAVE" } 


Check_deep = ["CM_argentina",
"CM_china",
"CM_colombia",
"CM_guinea",
"CM_guinea2",
"CM_caritro2",
"CM_NEUROBLASTOMA"]

Different_path = ["PernigoniN_2021", "MetaCardis_2020_a"]



Location = Path("/shares/CIBIO-Storage/CM/scratch/data/meta/")
Location2 = Path("/shares/CIBIO-Storage/CM/cmstore/data/meta/")

Missing_cohorts = []
Sample_files = {}
Missing = { "MP4":[], "SP4":[] }
for Cohort in  Dic_sample_project:
    print(Cohort)
    if Cohort in Rename: Cohort_work = Rename[Cohort]
    elif Cohort in Check_deep: Cohort_work = "Valles-ColomerM_2023_a/"+Cohort
    else: Cohort_work = Cohort
    
    if Cohort in Different_path:
        Location_cohort = Location2 / Path(Cohort_work)
    else:
        Location_cohort = Location / Path(Cohort_work)
    if not Location_cohort.exists(): Missing_cohorts.append(Cohort) ; continue
    
    MP4 = Location_cohort / "metaphlan-4.beta.1_vJan21_CHOCOPhlAnSGB_202103"
    SP4 = Location_cohort / "strainphlan-4beta_vJan21_CHOCOPhlAnSGB_202103"
    
    for Sample in Dic_sample_project[Cohort]:
            S_mp4 = MP4 / Sample
            S_sp4 = SP4 / Sample
            
            MP4_profile = S_mp4 / (Sample + "_unclassified.tsv")
            SP4_profile = S_sp4 / (Sample + ".pkl")

            if not MP4_profile.exists(): Missing["MP4"].append( (Sample, Cohort) )
            if not SP4_profile.exists(): Missing["SP4"].append( (Sample, Cohort) )

            Sample_files[Sample] = [Cohort, str(MP4_profile),str(SP4_profile)] 

with open("Missing.txt", "w") as I:
    Write = ""
    for i in Missing["MP4"]:
        Write += i[0] +"\t" + i[1] + "\tMetaphlan" + "\n"
    for i in Missing["SP4"]:
        Write += i[0] +"\t" + i[1] + "\tStrainphlan" + "\n"
    I.write(Write)

with open("File_location.txt", "w") as I:
    Write = ""
    for Sample in Sample_files:
            Info = Sample_files[Sample]
            Write += str(Sample) + "\t" + "\t".join(Info) + "\n"
    I.write(Write)    
    


