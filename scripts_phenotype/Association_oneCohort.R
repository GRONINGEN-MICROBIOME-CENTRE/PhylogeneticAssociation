library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(ape)
library(cmdstanr)
#set_cmdstan_path("/groups/umcg-gastrocol/tmp01/Sergio/.cmdstan/cmdstan-2.32.2/")
library(anpan)
library(furrr)
library(phytools)
library(tidyverse)

plan(multisession, workers = 4)


#Disease = "CRC"
#SGB = "t__SGB4933_group"
#Tree="/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Symlink_phylo/RAxML_bestTree.t__SGB4269.TreeShrink.tre"

source("Functions.R")


print("============== Initiating association script =============")
print("Usage: Rscript Association.R SGB_NAME PHYLOGENETIC_TREE [PHENOTYPE] [COV1,COV2,COV3...] [Location] [Cohort] ")
print("Checking variables included")

args <- commandArgs(trailingOnly = TRUE)
print(args)
if (!length(args) > 1) {
  message("Insufficient command-line arguments.")
  q() 
}

for (n in seq(length(args)) ){
  arg <- args[n]
  if (arg == "NA") {
    args[n] = NA
 }
}

SGB = args[1]
Tree = args[2]
Disease = args[3]
D = args[5] #D = "Assocaiton_per_cohort"
Covariates = args[4]
Covariates = str_split(Covariates, ",")[[1]] 
Cohort = args[6]
if (grepl(",",Cohort) ){
Cohort = str_split(Cohort, ",")[[1]]
} else {
Cohort = c(Cohort)
}

Do_Remove_repeated = T
Metadata_location ="../Phenotypes/Phenotypes_merged_complete.tsv" #"Phenotypes_merged.tsv"

print( paste0("Analyzing SGB: ", SGB, " and its tree: ", Tree, " in cohort: ", Cohort ) ) 

#1. Read Tree
print("Reading tree")
read.tree(Tree) -> Tree
if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }
Tree = midpoint.root(Tree)
#2. Rename labels tree
print("Renaming labels")
Clean_names =  function(Name){
  'Samples from 500FG FSK need everything after _ to be removed. If _metaphlan4 is attached at the end of the name, remove'
  if ( grepl("HV", Name)) {
    Name = str_split(Name, "_")[[1]][1]
  }   
  Name = str_replace(Name, "_metaphlan4", "") 
  return(Name)
}
Tree$tip.label %>% sapply(. , Clean_names) -> New_names
Tree$tip.label = New_names

#3. Get metadata
print("Reading metadata")
#Remove some datasets/. Sample MH0012 is repeated with the same name in two datasets. Study ThomasAM_2019_c and 500FG_FSK are repeated.
Metadata = read_tsv(Metadata_location) %>% filter(! study_name %in% c("ThomasAM_2019_c","500FG_FSK")  )
Metadata %>% group_by(ID_anal) %>% sample_n(1) -> Metadata
#Trimming metadata
Metadata %>% filter(ID_anal %in% Tree$tip.label) -> Metadata

Metadata %>% filter(study_name %in% Cohort)  -> Metadata
Cohort = paste(Cohort, collapse=",")


#4. Pruning tree
print("Pruning tree")
pruned_tree <- keep.tip(Tree, Metadata$ID_anal )
print(paste0("pruned_tree tips: ", length(pruned_tree$tip.label)))
#5. Sort Metadata
Metadata %>% filter(ID_anal %in% pruned_tree$tip.label) -> Metadata
Metadata %>% group_by(ID_anal) %>% summarise(N = n()) %>% arrange(desc(N))
Metadata <- arrange(Metadata, match(ID_anal, pruned_tree$tip.label ))
print(paste0("Metadata matched ", dim(Metadata)[1] ))
print(Metadata %>% group_by(study_name) %>% summarise(N = n()) %>% arrange(desc(N))  )
#6. Removal of repeated samples
#Usually they cluster together. We could in principle leverage some random effects to accomodate multiple 
if (Do_Remove_repeated == T){
	Unrepeated = Remove_repeated(Metadata, pruned_tree)
	Metadata = Unrepeated[[1]] ; pruned_tree = Unrepeated[[2]]
}


#8. Fit model
set.seed(890)

#Plot_name = paste0("Results/General_plots/Posterior_and_covariates/Disease/", paste(Covariates, collapse=","),"/", paste0(SGB,"_",Disease) ) 


Fit_model = function(Meta = Metadata, Phenotype_name=Disease, Tree=pruned_tree, COV=Covariates, Dir = NA ){
Results = tibble()
#Fam = "binomial"
#if(Phenotype_name %in% c("LDL", "HDL","HbA1c","triglycerides", "bilirubin",'SBP','DBP','HeartRate', "Age", "Sex", "BMI" ) ){ Fam = "gaussian" }
 
P = Meta[Phenotype_name]
P_noNA = P[!is.na(P)]
if (all(P_noNA %in% c(0, 1))) {
    Fam = "binomial"
} else { Fam = "gaussian" }


   
if (!file.exists( paste0("Results/",Dir))) { dir.create(paste0("Results/", Dir)) }
if (!file.exists( paste0("Results/",Dir,"Models" ))) { dir.create(paste0("Results/", Dir, "Models" )) }
if (!file.exists( paste0("Results/",Dir,"Plots" ))) { dir.create(paste0("Results/", Dir, "Plots" )) }
if (!file.exists( paste0("Results/",Dir,"Reports" ))) { dir.create(paste0("Results/", Dir, "Reports" )) }        
    
    
Plot_dir = paste0("Results/",Dir,"Plots/",SGB)
if (!file.exists(Plot_dir)) {
  dir.create(Plot_dir)
}
Summary_stats = paste0( "Results/", Dir, "Reports/",SGB, ".tsv")
Model_dir = paste0("Results/", Dir, "Models/",SGB)
if (!file.exists(Model_dir)) {
  dir.create(Model_dir)
}
Model_save = paste0(Model_dir, "/", Phenotype_name, "/" )
#if ( file.exists(paste0(Model_save, "Model.rds") ) ){ return() }


if (class(as_vector(Meta[Phenotype_name])) != "logical" & Fam=="binomial" ){
	Meta[Phenotype_name] =  ifelse(as.factor(as_vector(Meta[Phenotype_name])) == "1", T, F)
}
Labels = Meta$ID_anal[!is.na(Meta[Phenotype_name] ) ]
Meta %>% filter(ID_anal %in% Labels) -> Meta
keep.tip(Tree, Meta$ID_anal) -> Tree

Meta %>% mutate(sample_id = ID_anal, .before=1) %>% ungroup() %>% select("sample_id", "study_name", COV, Disease) %>% drop_na() -> Meta

#Check if multiple countries are avilable 
COV2 = c()
for (Covariate in COV){
        Meta[Covariate] %>% unique() %>% dim() -> Dim
        if (Dim[1] > 1 ){ COV2 = c(COV2, Covariate) }
}
COV2 = COV2[!COV2==Phenotype_name]
if (length(COV2) == 0){ COV2 = NULL }


anpan_pglmm(Meta, Tree, outcome = Phenotype_name, covariates = COV2, omit_na = T, family=Fam, out_dir=Model_save, save_object =T,parallel_chains = 4, show_plot_cor_mat=F, show_plot_tree=F,show_post=F) -> result
write_rds(result, paste0(Model_save, "Model.rds"))
print("Saving pareto")
#Save pareto
Do_pareto(result, paste0(Model_save, "pareto.tsv"))
print("Saving phylogenetic estimate")
Get_phylo_estimate(result, paste0(Model_save, "phylo_estimates.tsv"))

as_tibble(result$loo$comparison)[2,] %>% mutate(Phenotype = Phenotype_name, covariates= paste(COV2, collapse="+")  ) -> Result
Summary_stats_temp = paste0("Results/",Dir,"Reports/",SGB,"_", Phenotype_name,  ".tsv")
Result %>% mutate(SGB = SGB) -> Result
Result  %>% write_tsv(Summary_stats_temp)

Plot = paste0("Results/General_plots/Posterior_and_covariates/", Dir ,"/", SGB,"_",Phenotype_name) 
directory_path = paste0("Results/General_plots/Posterior_and_covariates/", Dir ,"/")
if (!dir.exists(directory_path)){
	dir.create(directory_path, recursive = TRUE)
}
Tree_plot_v2(SGB, Pheno=Phenotype_name, Plot_name=Plot,Model_name=result, Cov=COV, Tree=Tree)


}



D2 = paste0(D, "/",paste(Covariates, collapse=","),"_",Cohort, "/" )
M = Metadata
Tr = keep.tip(pruned_tree, M$ID_anal)		
C = Covariates
Fit_model(Meta = M, Phenotype_name = Disease, Tree= Tr, COV= C, Dir = D2)





