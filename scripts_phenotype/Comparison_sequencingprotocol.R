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

plan(multisession, workers = 4)


source("Functions.R")

args <- commandArgs(trailingOnly = TRUE)
print(args)
if (!length(args) == 2) {
  message("Insufficient command-line arguments.")
  q() 
}

SGB = args[1]
Tree = args[2]
Metadata_location = "../Phenotypes/Phenotypes_merged_complete.tsv"

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
#Remove some datasets/. Sample MH0012 is repeated with the same name in two datasets. Study ThomasAM_2019_c is repeated.
Metadata = read_tsv(Metadata_location) %>% filter(! study_name == "ThomasAM_2019_c" )
#Metadata %>% group_by(ID_anal) %>% sample_n(1) -> Metadata
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

#6. Prepare for comparison between FSK and APK 
Metadata %>% filter(study_name %in% c("500FG_FSK", "SchirmerM_2016") ) -> Metadata


Metadata %>% group_by(subject_id) %>% summarise( N = n() ) %>% arrange(desc(N))  %>% filter(N==2) -> Repeated_samples

Metadata %>% filter(subject_id %in% Repeated_samples$subject_id ) -> Metadata
pruned_tree <- keep.tip( pruned_tree, Metadata$ID_anal )
#Metadata %>% select( study_name, sample_id, id_used, subject_id, ID_anal) %>% print()

if (dim(Metadata)[1] < 20){ print("Few samples to compare tree") ; print(dim(Metadata)[1]) ; q() }

print(pruned_tree$tip.label)
order_indices <- match(Metadata$ID_anal, pruned_tree$tip.label)
Metadata = Metadata[order_indices, ]

Metadata %>% mutate(Protocol = ifelse(study_name == "500FG_FSK", 1, 0),ID_anal_p = ID_anal ,ID_anal=ifelse(study_name == "500FG_FSK", paste0(subject_id,"_FSK"),  paste0(subject_id,"_APK")   ) ) -> Metadata
pruned_tree$tip.label = Metadata$ID_anal

Metadata %>% group_by(Protocol) %>% summarise(n())

Subfolder = "Comparison_SequencingProtocol"
Summary_stats_temp = paste0("Results/",Subfolder ,"/Reports/",SGB, ".tsv")
Plot_d = paste0("Results/",Subfolder,"/Plots/", SGB, ".pdf" )
Model_dir = paste0("Results/", Subfolder, "/Models/",SGB)
if (!file.exists(Model_dir)) {
  dir.create(Model_dir)
}
Save_model = Model_dir

#8. Fit model
set.seed(890)
Results = tibble()

if (file.exists(Summary_stats_temp)) { q() }                   
                        
print(Save_model);print(Summary_stats_temp);print(Plot_d)
Result_tem = tryCatch({
             Function_fit_anpan("Protocol", pruned_tree, Metadata, Fam="binomial", Model_save = Save_model, Plot= Plot_d )  }, error = function(e){
             	Message = tibble( Phenotype = "Protocol",SGB =  SGB, Error = e$message)
             	Error = paste0("Results/",Subfolder ,"/Errors/",SGB,  ".tsv")
             	if (! file.exists(paste0("Results/",Subfolder ,"/Errors/"))){ dir.create(paste0("Results/",Subfolder ,"/Errors/") ) }
            	 print("Error!")
             	print(Message)
            	 write_tsv(Message, Error)
             	return(NULL)
             	})
             if ( is.null(Result_tem) ){ q() } 
             
Result_tem %>% mutate(SGB = SGB) -> Result_tem
Result_tem  %>% write_tsv(Summary_stats_temp)
#Tree_plot_v2(SGB, Pheno="Protocol")









