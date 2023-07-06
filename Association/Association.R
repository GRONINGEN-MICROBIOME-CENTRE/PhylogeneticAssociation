library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(ape)
library(cmdstanr)
set_cmdstan_path("/groups/umcg-gastrocol/tmp01/Sergio/.cmdstan/cmdstan-2.32.2/")
library(anpan)
library(furrr)

plan(multisession, workers = 4)


source("Functions.R")


print("============== Initiating association script =============")
print("Usage: Rscript Association.R SGB_NAME PHYLOGENETIC_TREE")
print("Checking variables included")

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) == 2) {
  message("Insufficient command-line arguments.")
  q() 
}

SGB = args[1]
Tree = args[2]
Do_Remove_repeated = T
Metadata_location = "../Phenotypes/Phenotypes_merged.tsv"
Variables_to_remove = "../Phenotypes/Remove_variables.tsv"

print( paste0("Analyzing SGB: ", SGB, " and its tree: ", Tree ) ) 

#1. Read Tree
print("Reading tree")
read.tree(Tree) -> Tree
if (class(Tree) == "multiPhylo" ){ Tree = Tree[[1]] }
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
#Remove some datasets/. Sample MH0012 is repeated with the same name in two datasets. Study ThomasAM_2019_c is repeated.
Metadata = read_tsv(Metadata_location) %>% filter(! study_name == "ThomasAM_2019_c" )
Metadata %>% group_by(ID_anal) %>% sample_n(1) -> Metadataread_tsv(Variables_to_remove) -> Not_associate
colnames(select(Metadata, -one_of( c("ID_anal", Not_associate$Remove)))) -> For_association
#Trimming metadata
Metadata %>% filter(ID_anal %in% Tree$tip.label) -> Metadata

#4. Pruning tree
print("Pruning tree")
pruned_tree <- keep.tip(Tree, Metadata$ID_anal )
print(paste0("pruned_tree tips: ", length(pruned_tree$tip.label)))
#5. Sort Metadata
Metadata %>% filter(ID_anal %in% pruned_tree$tip.label) -> Metadata
Metadata %>% group_by(ID_anal) %>% summarise(N = n()) %>% arrange(desc(N))
Metadata <- arrange(Metadata, match(ID_anal, pruned_tree$tip.label ))
print(paste0("Metadata matched ", dim(Metadata)[1] ))

#6. Removal of repeated samples
#Usually they cluster together. We could in principle leverage some random effects to accomodate multiple 
if (Do_Remove_repeated == T){
	Unrepeated = Remove_repeated(Metadata, pruned_tree)
	Metadata = Unrepeated[[1]] ; pruned_tree = Unrepeated[[2]]
}
#7. Name output
Plot_dir = paste0("Results/Plots/",SGB)
if (!file.exists(Plot_dir)) {
  dir.create(Plot_dir)
}
Summary_stats = paste0("Results/Reports/",SGB, ".tsv")
Model_dir = paste0("Results/Models/",SGB)
if (!file.exists(Model_dir)) {
  dir.create(Model_dir)
}
#8. Fit model
set.seed(890)
Results = tibble()

#Iteration through phenotypes and selection based on:
#Number of samples without NA >= 20
#If binary, make sure the labels are binary. If the number of cases is below 15, not run
Threshold_no_NA = 20
Threshold_n_cases = 15
Fit_anpan = T

print("Iterating through phenotypes")
for ( Phenotype in For_association ){
  #print(Phenotype)
  P = Metadata[Phenotype]
  P_noNA = P[!is.na(P)]
  if (length(P_noNA) < Threshold_no_NA){ next }
  if (all(P_noNA %in% c(0, 1))) {
    Family = "binomial"
    Counts = table(P)
    if (length(Counts) < 2) { next }
    if (Counts["1"] < Threshold_n_cases){ next }
  } else {
    Family = "gaussian"
  }
  print(paste0(Phenotype, " for association") )
  Save_model = paste0(Model_dir, "/", Phenotype, "/" )
  if (!file.exists(Save_model)) {
    dir.create(Save_model)
  }
  if (Fit_anpan == T){
	Plot_path =  paste0(Plot_dir, "/", Phenotype, ".pdf" )
	if (!file.exists(Plot_path)) { next }
  	Function_fit_anpan(Phenotype, pruned_tree, Metadata, Fam=Family, Model_save = Save_model, Plot= Plot_path  ) %>% rbind(Results, . ) -> Results
  } else {
	if (Remove_repeated == T){
		Function_fit_phyr(Phenotype, pruned_tree, Metadata, Fam=Family, Plot= paste0(Plot_dir, "/", Phenotype, "_ggtree.pdf" )  ) %>% rbind(Results, . ) -> Results
  	} else { 
		Function_fit_phyr(Phenotype, pruned_tree, Metadata, Fam=Family, Plot= paste0(Plot_dir, "/", Phenotype, "_ggtree.pdf" ), ID_random = "(1|subject_id)"  ) %>% rbind(Results, . ) -> Results
	}	
  }
}
Results %>% mutate(SGB = SGB) %>% write_tsv(Summary_stats)



