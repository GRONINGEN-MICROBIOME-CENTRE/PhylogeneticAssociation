library(tidyverse)
library(ape)
library(vegan)
library(phytools)

#Running with apptainer:
#apptainer shell --writable --no-home  --bind /scratch/p300317/Strains_sergio/:/mnt/ rstudio_instance/
#cd /mnt/Analysis/2_DistanceComparison
#Rscript Match_and_correlate.R

Prepare_table = function(){

    #1. Get Trees for analysis
    Sizes = read_tsv('Tree_sizes2.txt')
    Sizes %>% filter(N_tips >= 300) -> Sizes
    Sizes %>% mutate(Location_tree = paste0("Phylogenies/", Tree) ) -> Sizes
    #2. For each species, get the GTPro file, if exists
    ##There are SGB groups, in cases where we have a group, match to any of the SGBs within the group
    Groups = read_delim('../SGB_groups.csv', ';', col_types = "cc" )
    Groups_list = list()

    for (Entry in Groups$main_sgb) {
        Values = filter(Groups, main_sgb == Entry)$merged_sgbs
        Entry = paste0("t__", Entry)
        
        if (grepl(",", Values)) {
            Values_split = str_split(Values, ",")[[1]]
            Values_split = sapply(Values_split, function(x) paste0("t__SGB", x))
            Groups_list[[Entry]] = Values_split
            for (Name in Values_split){
                Groups_list[[Name]] = c(Entry, Values_split )
            } 
        } else {
            Name = paste0("t__SGB", Values)
            Groups_list[[Entry]] = c(Name)
            Groups_list[[Name]] = c(Entry)
        }
    }

    Location_GTPro = "../1_AssignTaxonomy/Assignment/Genomes_assigned/" # t__SGB4828.RData
    df_distances = tibble()

    for (Entry in Sizes$Tree) {
        SGB = str_split(Entry, "\\.")[[1]][2]
        
        if (grepl('group', SGB)) {
            SGB2 = str_replace(SGB, "_group", "")
            Entries = Groups_list[[SGB2]]
            Location_distance = NA
            for (Entry in Entries) {
                Location_distance2 = paste0(Location_GTPro, Entry, ".RData")
                Exists = file.exists(Location_distance2)
                if (Exists == FALSE) {
                    next
                } else { Location_distance =  Location_distance2}
            }
        } else {
            Location_distance = paste0(Location_GTPro, SGB, ".RData")
            Exists = file.exists(Location_distance)
            if (Exists == FALSE) {
                Location_distance = NA
            }
        }
        
        df_distances = df_distances %>% rbind(tibble(SGB = SGB, Location_distance = Location_distance))
    }

    Sizes$Tree  %>% sapply( function(Tree){ str_split(Tree, "\\.")[[1]][2] } ) -> Sizes$SGB
    left_join(Sizes,df_distances, by='SGB') -> SGB_info
    print(SGB_info)
    SGB_info %>% filter(is.na(Location_distance)) %>% print()
    SGB_info %>% drop_na() -> SGB_info
    write_tsv(SGB_info, 'InfoLocations.tsv')
}


##Useful fucntions for processing
Process_tree = function(Tree){
  read.tree(Tree) -> Tree
  if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }
  Tree = midpoint.root(Tree)
  #2. Rename labels tree
  print("Processing tree")
  Tree$tip.label %>% sapply(. , Clean_names) -> New_names
  Tree$tip.label = as.vector(New_names)
  if ( length(Tree$tip.label) < 100 ) { print("Tree is too small") ; q() }
  return(Tree)
  
}
Clean_names =  function(Name){
  'Samples from 500FG FSK need everything after _ to be removed. If _metaphlan4 is attached at the end of the name, remove'
  if ( grepl("HV", Name)) {
    Name = str_split(Name, "_")[[1]][1]
  }   
  Name = str_replace(Name, "_metaphlan4", "") 
  return(Name)
}
Get_distance = function(Location){
    load(Location)
    return(snp.dist)
}
Do_Test = function(SGB_info, Species){
    print(paste0("Running species ", Species) )
    Entry = filter(SGB_info, SGB == Species )
    print('Reading and processing phylogenetic tree')
    Tree = Process_tree(Entry$Location_tree)
    print('Reading SNP distance') 
    SNP_Distance = Get_distance(Entry$Location_distance) 
    #Get distance from Tree
    Phylo_dist = cophenetic.phylo(Tree) %>% as.matrix()
    #Filter samples
    print('Matching data')
    intersect(colnames(Phylo_dist), colnames(SNP_Distance)) -> Samples
    Phylo_dist = Phylo_dist[Samples, Samples]
    SNP_Distance = SNP_Distance[Samples, Samples]
    #Mantel
    print('Running Mantel')
    mantel( as.dist(Phylo_dist), as.dist(SNP_Distance), method = "pearson", permutations=1000 ) -> Test
    N_samples = length(Samples)
    tibble(SGB=Species,Rho= Test$statistic, P=Test$signif, Permutations= Test$permutations, N= N_samples ) -> Results
    return(Results)
}


#3. Read tree and make into a distance, read GTPro distances
Run_mantel = function(SGB_info){
    #Tree    N_tips Location_tree SGB   Location_distance
    Result_all = tibble()
    #4. Run mantel test
    for (Species in SGB_info$SGB){
        Results = Do_Test(SGB_info, Species)
        rbind(Results_all, Results) -> Results_all
    }
    write_tsv(Results_all, 'Correlations_per_taxa.tsv')
}

Run_mantel_paral = function(SGB_info, SGB){
    Out = paste0("Results_mantel/", SGB, ".tsv")
    Results = Do_Test(SGB_info, SGB)
    write_tsv(Results, Out)
}

#Generate info location
if (file.exists( 'InfoLocations.tsv') ){
    SGB_info = read_tsv( 'InfoLocations.tsv')
} else {
    SGB_info = Prepare_table()
}


#Receive SGB name from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    print("No arguments given, running loop for all SGB")
    Run_mantel(SGB_info)
} else {
    S = args[1]
    print( paste0("Running with SGB: ", S ) )
    Run_mantel_paral(SGB_info, S)
}




