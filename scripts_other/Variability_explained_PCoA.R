###PCA####
library(tidyverse)
source("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Make_Associations/Association/Functions.R")
library(vegan)
library(phytools)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2){
	print("No command line arguments given. Using default SGB and tree")
	SGB = "SGB4837_group"
	Tree = "../../Symlink_phylo/RAxML_bestTree.t__SGB4837_group.TreeShrink.tre"
}
SGB = args[1]
Tree = args[2]
Metadata_location ="../Phenotypes/Phenotypes_merged_complete.tsv" #"Phenotypes_merged.tsv"
Only_healthy = T


if (Only_healthy == T){
	Output = paste0("Results/PCA_partion/", SGB, "_noDisease_onlyHealthy.tsv")
} else { Output = paste0("Results/PCA_partion/", SGB, "_noDisease.tsv") }

if (file.exists(Output)){
	print("Output file exists, exiting script")
}

print("Processing tree")

#1. Read. Tree
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

#3. Read metadata
Metadata = read_tsv(Metadata_location) %>% filter(! study_name %in% c("ThomasAM_2019_c","500FG_FSK")  )
Metadata %>% group_by(ID_anal) %>% sample_n(1) -> Metadata
Metadata %>% filter(ID_anal %in% Tree$tip.label) -> Metadata

if (Only_healthy == T){
	Metadata %>% filter(disease == "healthy" | ( is.na(disease) & NoDisease==1  ) ) -> Metadata
}

pruned_tree <- keep.tip(Tree, Metadata$ID_anal )

#3.Removal of repeated samples
#Usually they cluster together. We could in principle leverage some random effects to accomodate multiple 
Unrepeated = Remove_repeated(Metadata, pruned_tree)
Metadata = Unrepeated[[1]] %>% ungroup() ; pruned_tree = Unrepeated[[2]]


print("Making a pcoa (mds) from phylogenetic distance matrix")

#4. Make PCA
cophenetic.phylo( pruned_tree ) -> Phylo_dist
pcoa = cmdscale (Phylo_dist, eig = TRUE, k = dim(Metadata)[1]-1  )

Eigen_values = pcoa$eig / (pcoa$eig %>% sum())
Cumulative_variability_explained = cumsum(Eigen_values) 
PCs_check = which.max(Cumulative_variability_explained >= 0.9)

print( paste0("Number of PCs to reach 90% variability explained: ", PCs_check) )


pcoa$points %>% as.data.frame() %>% rownames_to_column("ID_anal") %>% as_tibble() -> Components

Components %>% left_join(Metadata) %>% filter(! is.na(Age) )  -> Components


print("Getting list of countries with more than one study")
#Find countries with more than 1 study
Components %>% group_by(study_name, Country) %>% summarise(n()) %>% group_by(Country) %>% summarise(N = n())  %>% filter(N > 1) -> countries_multiple
#Find continents with more than 1 study
Components %>% group_by(study_name, Continent) %>% summarise(n()) %>% group_by(Continent) %>% summarise(N = n())  %>% filter(N > 1) -> continents_multiple
#Remove samples from studies with few samples
Components %>% group_by(study_name) %>% summarise(N = n()) %>%  filter(N < 10) -> studies_remove

print(continents_multiple)
Components %>% filter(! study_name %in% studies_remove$study_name ) %>% filter(Country %in% countries_multiple$Country ) %>% filter(! (is.na(Age) | is.na(BMI) | is.na(Continent) | is.na(Sex) ) ) %>% group_by(Country, Continent, study_name) %>% summarise(n()) -> Final



if (dim(Final)[1] < 3 | length(unique(Final$study_name)) < 2  ){
	file.create(Output)
	q()
}	


All_results = tibble()



print("Fitting GLMM with nested random effects, using MCMC")

Use_MCMC = function(){
library(MCMCglmm)

prior1 = list(
  G=list(
    G1=list(V=1,nu = 0.001),
    G2=list(V=1,nu=0.001),
    G3=list(V=1,nu=0.001)),
  R=list(V=1,nu=0.001))




Get_info_from_MCMC = function(VCV){
  as_tibble(VCV) %>% summarise_all(median) -> info_median
  as_tibble(VCV) %>% summarise_all(mean) -> info_mean
  colnames(info_median) = paste0("median_", colnames(info_median) )
  colnames(info_mean) = paste0("mean_", colnames(info_mean) )
  posterior.mode(VCV) %>% t() %>% as.data.frame() %>% as_tibble() -> info1
  HPDinterval(VCV,prob = 0.95) %>% as.data.frame() %>% rownames_to_column('Variable')  %>%
    pivot_wider(names_from = Variable, values_from = c(lower, upper)) -> info2
  cbind(info1, info2) %>% cbind(info_median) %>% cbind(info_mean) %>% as_tibble()  %>% return()
}

All_results = tibble()
Components %>% filter(! study_name %in% studies_remove$study_name ) %>% filter(Country %in% countries_multiple$Country ) %>%
  select( - family ) %>% filter(!is.na(Age) ) %>% filter(!is.na(Sex) ) %>% filter(!is.na(BMI) ) -> Input
print(Input)
for (N in seq(PCs_check)){
  print(paste0("Running component ", N) )
  Formula = as.formula(paste0( "V",N, " ~ Age + Sex + BMI"))
  Input %>%  MCMCglmm(Formula ,
           random = ~Continent + Continent:Country+ Continent:Country:study_name, 
           prior = prior1,
           thin=50,
           data=. ,
           family = "gaussian",nitt = 60000,burnin=10000 ) -> result
  VCV2 = result$VCV
  VCV2 / rowSums(VCV2) -> VCV3
  #plot(VCV3)
  
  Get_info_from_MCMC(VCV2) %>% mutate(Component = N, .before=1) %>%  rbind(All_results, . ) -> All_results
}

write_tsv(All_results %>% mutate(N_samples = dim(Input)[1])  , Output )

}

Use_lmer = function(){
library(lmerTest)
Input = Components %>% filter(! study_name %in% studies_remove$study_name ) %>% filter(Country %in% countries_multiple$Country ) %>% filter(! (is.na(Age) | is.na(BMI) | is.na(Continent) | is.na(Sex) ) )
for (N in seq(PCs_check)){
  Component = paste0("V", N )
  
  Formula2 <- paste0(Component, " ~ Age + Sex + BMI")

# If multiple continents
if (length(unique(Final$Continent)) > 1) {
  Formula2 <- paste0(Formula2, " + (1|Continent)")
  # If multiple countries
  if (length(unique(Final$Country)) > 1) {
    Formula2 <- paste0(Formula2, " + (1|Continent:Country) + (1|Continent:Country:study_name)")
  }
}

# If only one continent
if (length(unique(Final$Continent)) == 1) {
  # If multiple countries
  if (length(unique(Final$Country)) > 1) {
    Formula2 <- paste0(Formula2, " + (1|Country:study_name)")
  }
}

# If multiple continents, but one country per continent
Final %>% group_by(Continent) %>% summarise(N = n()) -> Number_countries_per_continent
if (dim(Number_countries_per_continent)[1] > 1 & dim(Number_countries_per_continent)[1] == dim(filter(Number_countries_per_continent, N == 1))[1]) {
  Formula2 <- paste0(Formula2, " + (1|Continent:study_name)")
} 


  Input %>% lmer(Formula2, . ) -> Result
  
  tryCatch({
  CI <- confint(Result)
  }, error = function(e) {
  # If an error occurs, calculate confidence intervals using confint.merMod(model, method = "Wald")
  CI <- confint.merMod(Result, method = "boot")
})

  Result %>% summary() -> Result
  Result$varcor %>% as_tibble() %>% cbind(CI[1:dim(.)[1],]) %>% as_tibble() %>% select(-c(var1,var2)) %>% rename(Random_effect = grp) -> RE

  #Result$varcor %>% cbind(CI[1:3,] ) %>% as.data.frame() %>% rownames_to_column("Random_effect")  -> RE
  #colnames(RE)[2] = "Variance_Explained"
  RE %>% mutate(Component= N, .before=1) %>% rbind(All_results, . ) -> All_results
  #Result$coefficients %>% as.data.frame() %>% rownames_to_column("Features") %>% as_tibble() %>% mutate(Component= Component) %>% rbind(All_results, . ) -> All_results
}
left_join(All_results, data.frame(Component = seq(N), Variance_explained_PC = Eigen_values[1:N])) -> All_results
All_results[] <- lapply(All_results, unlist)

write.table(All_results %>% mutate(N_samples = dim(Input)[1])  , Output, quote=F, sep="\t", row.names =F )

}          


Use_lmer()


