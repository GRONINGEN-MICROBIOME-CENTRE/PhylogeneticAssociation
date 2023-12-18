library(tidyverse)

###########################################
####Geography assocaitons and taxonomy#####
###########################################

SGB_taxonomy = read_tsv("/mnt/project/Make_Associations/Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2", col_names = F ) %>% rename(SGB=X1, Taxonomy=X2)
Geography_assocations = read_tsv("/mnt/project/Make_Associations/Association/Results/Geography_cor/Association_results_merged.tsv")
Geography_assocations %>% mutate(SGB = str_replace(SGB, "t__", "") ) -> Geography_assocations
SGB_taxonomy <- SGB_taxonomy %>%
  separate(Taxonomy, into = c("domain","phylum", "class", "order", "family", "genera", "species"), sep = "\\|", remove = FALSE)

left_join(Geography_assocations, SGB_taxonomy) ->  Geography_assocations
Geography_assocations %>% filter(FDR<0.05) %>%  gather(., Taxonomic_level, Taxa, domain:species, factor_key=TRUE) %>% filter(!Taxonomic_level %in% c("domain", "species", "class", "oder", "family") ) %>% 
  ggplot(aes(x=Taxonomic_level, fill=Taxa)) + geom_bar() + theme_bw()


Geography_assocations %>% group_by(FDR<0.05, phylum) %>% summarise(N = n()) %>% ungroup() %>%
  group_by(phylum) %>%
  summarise(
    p_value = fisher.test(matrix(c(sum(`FDR < 0.05`), sum(`FDR < 0.05` == FALSE), sum(N - `FDR < 0.05`), sum(`FDR < 0.05` == TRUE)), nrow = 2))$p.value
  ) -> Fisher_results

##Test differences of rho
library(fgsea)
taxa_list = list()
for (i in c("phylum", "class", "order", "family", "genera") ){
  Geography_assocations[i] %>% as_vector() -> Names
  taxa_list_pre <- split(Geography_assocations$SGB, Geography_assocations[i] %>% as_vector() %>% as.vector() )
  taxa_list <- c(taxa_list, taxa_list_pre)
}
Geography_assocations %>% arrange(desc(Rho))-> R1
#Prepare ranks
ranks = R1$Rho
names(ranks) = as.character(R1$SGB)
#Run GSEA
fgseaRes <- fgsea::fgsea(taxa_list, ranks)
as_tibble(fgseaRes) %>% arrange(padj) -> Assoc_rho_taxa
#plots
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
fgsea::plotGseaTable(taxa_list[topPathways], ranks, fgseaRes,  gseaParam=0.5)

plotEnrichment(taxa_list[["f__Lachnospiraceae"]], ranks) + labs(title="Lachnospiraceae") + theme( text = element_text(size = 21),  plot.title = element_text(size = 26), axis.title.x = element_text(size = 21),  axis.title.y = element_text(size = 21)  )
plotEnrichment(taxa_list[["o__Clostridiales"]], ranks) + labs(title="Clostridiales")

New_Edge = c()
for ( edge in  Assoc_rho_taxa$leadingEdge){
  New_Edge = c(New_Edge, paste(edge, collapse=","))
}
Assoc_rho_taxa %>% mutate(leadingEdge = New_Edge) %>% write_tsv("/mnt/project/Make_Associations/Association/Results/Geography_enrichment.tsv")


results_list <- list()

# Loop through unique phyla
unique_phyla <- unique(data$phylum)
for (phylum in unique_phyla) {
  # Subset data for the current phylum
  subset_data <- data[data$phylum == phylum, ]
  
  # Perform Fisher's exact test
  Matrix = matrix(c(sum(subset_data$`FDR < 0.05`), sum(!subset_data$`FDR < 0.05`), sum(subset_data$N - subset_data$`FDR < 0.05`), sum(subset_data$`FDR < 0.05`)), nrow=2)
  fisher_result <- fisher.test(, nrow = 2))
  
  # Store results in a list
  results_list[[phylum]] <- list(
    phylum = phylum,
    p_value = fisher_result$p.value
  )
}  


Run_enrichment = function(Geography_assocations, Level){
  # Subset data with FDR < 0.05
  subset_data <- Geography_assocations[Geography_assocations$FDR < 0.05, ]
  # Initialize an empty list to store results
  results_enrichement = tibble()
  # Loop through unique phyla in the subset
  unique_phyla <- unique(subset_data[Level] %>% as_vector() )
  for (phylum in unique_phyla) {
    # Create a contingency table
    #Column specific taxa
    N_phyla = sum(Geography_assocations[Level] %>% as_vector() == phylum)
    Sig = sum(subset_data[Level] == phylum)
    NoSig = N_phyla - Sig
    #Column all taxa
    N_phyla_other = sum(Geography_assocations[Level]  %>% as_vector() != phylum)
    Sig_other = sum(subset_data[Level]  %>% as_vector() != phylum)
    NoSig_other = N_phyla_other - Sig_other
    #Make matrix
    contingency_table <- matrix( c(Sig, NoSig, Sig_other,  NoSig_other), nrow = 2, ncol = 2 ) %>% as.table()

    #run test
    fisher_result <- fisher.test(contingency_table, alternative = "greater")
    #save
    results_enrichement %>% rbind(tibble(Taxa=phylum, P=fisher_result$p.value, N_taxa = N_phyla, Sig_taxa=Sig )) -> results_enrichement
  }
  
  results_enrichement %>% mutate(FDR = p.adjust(P, "fdr")) %>% return()
}

Results_enrichment = tibble()
for (i in c("phylum", "class", "order", "family", "genera") ){
  Run_enrichment(Geography_assocations, i) %>% mutate(Taxonomic_level = i ) %>% rbind(Results_enrichment, . ) ->  Results_enrichment
}
Results_enrichment %>% mutate(FDR_all = p.adjust(p, "fdr") )



###########################################
####Overall associations and taxonomy#####
###########################################
Overall_associations =  read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/All_results.csv" ) %>% mutate(ID = paste0(Phenotype,SGB) )
Significant_associations = read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/All_results_filtered.csv" )%>% mutate(ID = paste0(Phenotype,SGB) )
Overall_associations %>% mutate(FDR = ifelse(ID  %in% Significant_associations$ID, 0, 1 )) -> Overall_associations

Overall_associations %>%  mutate(SGB = str_replace(SGB, "t__", "") ) -> Overall_associations
left_join(Overall_associations, SGB_taxonomy) ->  Overall_associations

Results_enrichment_p = tibble()
for (i in c("phylum", "class", "order", "family", "genera",  "SGB") ){
  Run_enrichment(Overall_associations, i) %>% mutate(Taxonomic_level = i ) %>% rbind(Results_enrichment_p, . ) ->  Results_enrichment_p
}
Results_enrichment_p %>% mutate(FDR_all = p.adjust(P, "fdr") ) %>% arrange(P) -> Results_enrichment_p
write_tsv(Results_enrichment_p, "/mnt/project/Make_Associations/Association/Results/Summaries/Enrichment_taxonomiclevel_associations.tsv")

Results_enrichment_p %>% filter(grepl("SGB", Taxa)) %>% left_join(SGB_taxonomy %>% rename(Taxa = SGB) %>% select(Taxa, species) ) 



###########################################################################
####Enrichment of taxonomic level in continent-stratified analysis#########
##########################################################################
#############################################
####Summarise associations per continent#####
#############################################

read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Africa_filtered.csv") %>% mutate(Continent="Africa") %>% mutate(ID = paste0(Phenotype,"-",SGB) ) -> Africa_filtered
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Asia_filtered.csv") %>% mutate(Continent="Asia") %>% mutate(ID = paste0(Phenotype,"-",SGB) ) -> Asia_filtered
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Europe_filtered.csv") %>% mutate(Continent="Europe") %>% mutate(ID = paste0(Phenotype,"-",SGB) ) -> Europe_filtered
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica_filtered.csv") %>% mutate(Continent="North_America") %>% mutate(ID = paste0(Phenotype,"-",SGB) ) -> NorthAmerica_filtered
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_SouthAmerica_filtered.csv")  %>% mutate(Continent="South_America") %>% mutate(ID = paste0(Phenotype,"-",SGB) )  -> SouthAmerica_filtered
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Oceania_filtered.csv") %>% mutate(Continent="Oceania") %>% mutate(ID = paste0(Phenotype,"-",SGB) ) -> Oceania_filtered

#info phenotypes
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypic_categories.tsv") %>% rename(Phenotype = Phenos) -> Categories


All_associations = rbind(rbind(rbind(rbind(rbind(Africa_filtered, Asia_filtered), Europe_filtered), NorthAmerica_filtered), SouthAmerica_filtered), Oceania_filtered) 
All_associations %>% left_join(Categories) -> All_associations

#Global assocations
dim(All_associations)[1]
#Unique associations
All_associations %>% group_by(ID) %>% summarise(n()) %>% dim()
#Assocations in more than one continent
All_associations %>% group_by(ID) %>% summarise(N = n()) %>% filter(N > 1) %>% dim() 
#Number phenos
All_associations %>% group_by(Phenotype) %>% summarise(N = n())  %>% dim()
#Number SGB
All_associations %>% group_by(SGB) %>% summarise(N = n())  %>% dim()

#Per category (how many unique ID per category)
All_associations %>% group_by(Category) %>%  summarise(unique_ids = n_distinct(ID))
#To what
for (Cate in unique(All_associations$Category)){
  print(Cate)
  All_associations %>% filter(Category == Cate ) %>% group_by(Phenotype) %>% summarise(N = n_distinct(ID)) %>% arrange(desc(N)) %>% print()
  
}
#N Per continent
All_associations %>% group_by(Continent) %>% summarise(N = n()) %>% arrange(N)
#Are number of assoacitons related to number of samples?
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> Phenotypes
Phenotypes %>% group_by(Continent) %>% summarise(N_samples = n()) %>% left_join( All_associations%>% group_by(Continent) %>% summarise(N_assocations = n()), .  ) -> Comparison_Nvs_N
cor.test(Comparison_Nvs_N$N_assocations, Comparison_Nvs_N$N_samples, method = "spear")
#Are there taxonomic groups enriched?
All_associations %>%  distinct(ID, .keep_all = T) -> Significant_associations
Overall_associations =  read_csv("/mnt/project/Make_Associations/Association/All_results.csv" ) %>% mutate(ID = paste0(Phenotype,"-",SGB) )
Overall_associations %>% mutate(FDR = ifelse(ID  %in% Significant_associations$ID, 0, 1 )) -> Overall_associations
left_join(Overall_associations, SGB_taxonomy ) ->  Overall_associations
Results_enrichment_p = tibble()
for (i in c("phylum", "class", "order", "family", "genera",  "SGB") ){
  Run_enrichment(Overall_associations, i) %>% mutate(Taxonomic_level = i ) %>% rbind(Results_enrichment_p, . ) ->  Results_enrichment_p
}
Results_enrichment_p %>% mutate(FDR_all = p.adjust(P, "fdr") ) %>% arrange(P) -> Results_enrichment_p
write_tsv(Results_enrichment_p, "/mnt/project/Make_Associations/Association/Results/Summaries/Enrichment_taxonomiclevel_associations.tsv")
Results_enrichment_p %>% filter(grepl("SGB", Taxa)) %>% left_join(SGB_taxonomy %>% rename(Taxa = SGB) %>% select(Taxa, species) ) 
Results_enrichment_p %>% filter(FDR_all < 0.05) %>% filter(Taxonomic_level == "SGB") %>% rename(SGB = Taxa) %>% left_join(SGB_taxonomy) %>% select(X2)
Results_enrichment_p %>% filter(Taxonomic_level == "SGB")
#are those associating bcs sample size? or due to geographical confounding?
#sample size
read_rds("/mnt/project/SGB_abundances/Data/Prevalence.rds") -> Prevalence
Prevalence %>% select(-ID_anal) %>% apply(2, function(x){ sum(x)/length(x) } ) -> Prev
SGB_name = names(Prev) %>% lapply(function(x){str_split(x, "\\|")[[1]]-> y ; y[length(y)]  } ) %>% unlist()
tibble(Prevalence= Prev, SGB = SGB_name ) -> Prevalence
Results_enrichment_p %>% filter(Taxonomic_level == "SGB") %>% rename(SGB=Taxa) %>%  left_join(Prevalence) %>% lm(Sig_taxa ~ Prevalence, . ) %>% summary()
read_tsv("/mnt/project/Make_Associations/Association/Tree_sizes2.txt") -> Tree_size
Tree_size$Tree %>% sapply(function(x){ str_split(x, "\\.")[[1]][2] } ) -> SGB_names
Tree_size %>% mutate(SGB = SGB_names) -> Tree_size
Results_enrichment_p %>% filter(Taxonomic_level == "SGB") %>% rename(SGB=Taxa) %>%  left_join(Tree_size) -> Info_enrichment
Info_enrichment  %>% lm(Sig_taxa ~ N_tips, . ) %>% summary()
#geographical effect
read_tsv("/mnt/project/Make_Associations/Association/Results/Geography_cor/Association_results_merged.tsv") %>% left_join(Info_enrichment, . , by="SGB") %>% lm(Sig_taxa ~ N_tips + Rho, .) %>% summary()







