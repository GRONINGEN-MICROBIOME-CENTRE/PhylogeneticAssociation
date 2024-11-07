library(tidyverse)
library(patchwork)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")


Info = read_tsv("/mnt/project/Make_Associations/Association/Tree_sizes2.txt")
SGB_taxonomy = read_tsv("/mnt/project/Make_Associations/Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2", col_names = F )

#1. From Info, get the SGB name and retriet its complete taxonom from SGB_taxonomy.
extract_parts <- function(data_frame, column_name) {
  extracted_names <- sub(".*\\.(.*)\\.Tree.*", "\\1", data_frame[[column_name]])
  numeric_parts <- gsub("\\D+", "", extracted_names)
  # Extract part before numbers and after "__"
  before_numbers <- str_split(extracted_names, "_") %>% sapply(., function(x) x[[3]]) %>% gsub("\\d+", "", .)
  return(tibble(Tree =data_frame$Tree  ,SGB = str_replace(extracted_names, "t__","") , SGB_name = numeric_parts, Type = str_replace(before_numbers, "_group","")  ))
}
SGB_Tree = extract_parts(Info, "Tree")
left_join(Info, SGB_Tree) -> Info

#2. Assess whether the tree was used for analysis or not
Info %>% filter( Type == "SGB") %>% mutate(Included_analysis = ifelse(N_tips < 300, F, T ) ) -> Info
print("Check all trees used for association")
Info %>% group_by(Included_analysis ) %>% summarise(N = n()) %>% print()

#3. Combine tree info with taxonomy info
SGB_taxonomy %>% rename(SGB=X1, Taxonomy=X2 ) %>% mutate(Known = ifelse( grepl(Taxonomy, "_SGB"), F, T ) ) -> SGB_taxonomy
Info = left_join(Info,SGB_taxonomy)
#High level overview
Info %>% mutate(Phylum = str_split(Taxonomy, "\\|", simplify = TRUE)[, 2]) %>% mutate(Phylum = str_split(Phylum, "__", simplify = TRUE)[, 2])  -> Info
Info %>% rename(Number_samples_phylogeny = N_tips  ) -> Info

#4. Combine all info with SGB tree
read.tree("/mnt/project/Make_Associations/Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk") -> SGB_tree
keep.tip(SGB_tree, Info$SGB_name) -> SGB_tree

#5. Get information about: 
#       Does the tree show sequencing protocol differences?
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/SequencingProtocol_results.csv") -> Info_protocol
Info_protocol %>% select(SGB) %>% mutate( Show_protocol_differences = F ) %>%  mutate(SGB = str_remove(SGB, "t__")) -> Info_protocol
Info %>% left_join(Info_protocol) -> Info
#       Does the tree show geographical differences
read_tsv("/mnt/project/Make_Associations/Association/Results/Geography_cor/Association_results_merged.tsv") -> Info_geography
#read_csv("/mnt/project/Make_Associations/Association/Geography_results.csv") -> Info_geography
#read_csv("/mnt/project/Make_Associations/Association/Geography_results_filtered.csv") -> Info_geography_sign
#Info_geography %>% select(SGB) %>% distinct() %>% mutate(Geographical_effect = ifelse(SGB %in% Info_geography_sign$SGB, T, F )) %>%  mutate(SGB = str_remove(SGB, "t__")) -> Info_geography
Info_geography %>% mutate(Geographical_effect = ifelse(FDR<0.05, T, F ) ) %>% mutate(SGB = str_remove(SGB, "t__")) -> Info_geography

Info %>% left_join(Info_geography) -> Info
#       Is the tree associated to phenotypes? Quantify per phenotype group: e.g Antrophometric / Diseases / Biochemical / etc



#6. Make some plots

color_mapping <- c(
  "Actinobacteria" = c25[1],
  "Firmicutes" = c25[2],
  "Bacteroidetes" = c25[3],
  "Proteobacteria" = c25[4],
  "Euryarchaeota" = c25[5],
  "Synergistetes" = c25[7],
  "Lentisphaerae" = c25[8],
  "Kiritimatiellaeota" = c25[9],
  "Tenericutes" = c25[10],
  "Elusimicrobia" = c25[11],
  "Crenarchaeota" = c25[13],
  "Spirochaetes" = c25[14],
  "Fusobacteria" = c25[15],
  "Verrucomicrobia" = c25[15],
  "Planctomycetes" = c25[17],
  "Bacteria_unclassified" = c25[18],
  "Candidatus_Melainabacteria" = c25[19],
  "Candidatus_Thermoplasmatota" = c25[20],
  "Candidatus_Saccharibacteria" = c25[21]
)

Info %>%   select(SGB_name, everything())  %>% mutate(Log10_number_samples = log10(Number_samples_phylogeny)  ) %>% mutate(Taxonomic_group = ifelse(grepl("Archa", Taxonomy ), "Archaea", "Bacteria" ) ) -> Info2

#Reroot on archaea
Phylum_name = "Euryarchaeota"
leafs_ids =  filter(Info, Phylum == Phylum_name)$SGB_name
leaf_nodes <- match(leafs_ids, SGB_tree$tip.label)
Root = leaf_nodes[3] #2066
SGB_tree %>% root(., outgroup = Root) -> SGB_tree2

Find_node_ancestor = function(Phylum_name = "Firmicutes", subsample=50, SGB_tree=SGB_tree2){
  leafs_ids =  filter(Info, Phylum == Phylum_name)$SGB_name
  if (length(leafs_ids) > subsample){ sample(leafs_ids, subsample) -> leafs_ids } 
  leaf_nodes <- match(leafs_ids, SGB_tree$tip.label)
  #thresh = mean(leaf_nodes) + 3*sd(leaf_nodes)
  #leaf_nodes = leaf_nodes[ leaf_nodes < thresh   ]
  if (Phylum_name == "Proteobacteria"){ leaf_nodes = leaf_nodes[leaf_nodes<600] }
  if (Phylum_name == "Actinobacteria"){ leaf_nodes = leaf_nodes[leaf_nodes>1500] }
  
  common_ancestor <- getMRCA(SGB_tree, leaf_nodes)
  print(Phylum_name) ; print(leaf_nodes) ; print(common_ancestor)
  tibble(node =common_ancestor, Phylum_group = Phylum_name  ) %>% return()
}

set.seed(1343)
Find_node_ancestor("Actinobacteria", subsample=500) %>% rbind( Find_node_ancestor("Bacteroidetes") ) %>% 
  rbind( Find_node_ancestor("Firmicutes") ) %>% rbind(Find_node_ancestor("Proteobacteria",subsample=500 )) %>%
  rbind( Find_node_ancestor("Verrucomicrobia",subsample=200) ) %>% rbind(Find_node_ancestor("Lentisphaerae") ) %>% rbind( Find_node_ancestor("Candidatus_Melainabacteria") ) -> Annotation_phyla

Annotation_phyla

Info2

ggtree(SGB_tree2, layout="fan", open.angle=15, size=0.1)  %<+% (Info2 %>% mutate( Show_protocol_differences = ifelse(is.na(Show_protocol_differences), "NA", Show_protocol_differences ), Geographical_effect=ifelse(is.na(Geographical_effect), "NA", Geographical_effect ) ) )  -> p
p +   geom_balance(node=filter(Annotation_phyla, Phylum_group=="Bacteroidetes")$node  , fill=color_mapping["Bacteroidetes"], color=NA, alpha=0.3) +  #Bacteroidetes
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Firmicutes")$node  , fill=color_mapping["Firmicutes"], color=NA, alpha=0.3) + #Firmicutes
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Proteobacteria")$node  , fill=color_mapping["Proteobacteria"], color=NA, alpha=0.3) + #Proteobacteria
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Verrucomicrobia")$node  , fill=color_mapping["Verrucomicrobia"], color=NA, alpha=0.3) + #Verrucomicrobia
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Lentisphaerae")$node  , fill=color_mapping["Lentisphaerae"], color=NA, alpha=0.3) + #Lentisphaerae
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Candidatus_Melainabacteria")$node  , fill=color_mapping["Candidatus_Melainabacteria"], color=NA, alpha=0.3) +  #Candidatus_Melainabacteria
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Actinobacteria")$node  , fill=color_mapping["Actinobacteria"], color=NA, alpha=0.3) + #Actinobacteria
  #geom_tippoint(aes_string(col = "Phylum"  ), size= 0.7)  + scale_color_manual(values = color_mapping) + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Included_analysis), width=0.03,offset=0.1 ) + scale_fill_manual( values = c( "FALSE"="grey", "TRUE" = "#E83845") ) + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Log10_number_samples), width=0.03,offset=0.1 ) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Taxonomic_group), width=0.03,offset=0.1 ) + scale_fill_manual( values = c("steelblue", "darkgreen") ) + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Show_protocol_differences), width=0.03,offset=0.1 ) + scale_fill_manual( values = c("FALSE"="grey", "NA"="white" ) ) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Geographical_effect), width=0.03,offset=0.1 ) + scale_fill_manual( values = c("FALSE"="grey", "TRUE"="black", "NA"="white" ) ) -> SGB_Tree
  

#Version 2
ggtree(SGB_tree2, layout="fan", open.angle=15, size=0.1)  %<+% (Info2 %>% mutate( Show_protocol_differences = ifelse(is.na(Show_protocol_differences), "NA", Show_protocol_differences ), Geographical_effect=ifelse(is.na(Geographical_effect), "NA", Geographical_effect ) ) )  -> p
p +   geom_balance(node=filter(Annotation_phyla, Phylum_group=="Bacteroidetes")$node  , fill=color_mapping["Bacteroidetes"], color=NA, alpha=0.2) + geom_cladelab(node=filter(Annotation_phyla, Phylum_group=="Bacteroidetes")$node, label= "Bacteroidetes", angle=0, fontsize=3, vjust=.1  ) + #Bacteroidetes
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Firmicutes")$node  , fill=color_mapping["Firmicutes"], color=NA, alpha=0.2) + geom_cladelab(node=filter(Annotation_phyla, Phylum_group=="Firmicutes")$node, label= "Firmicutes", angle=0, fontsize=3, vjust=.1  ) + #Firmicutes
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Proteobacteria")$node  , fill=color_mapping["Proteobacteria"], color=NA, alpha=0.2) + geom_cladelab(node=filter(Annotation_phyla, Phylum_group=="Proteobacteria")$node, label= "Proteobacteria", angle=0, fontsize=3, vjust=.5  ) + #Proteobacteria
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Verrucomicrobia")$node  , fill=color_mapping["Verrucomicrobia"], color=NA, alpha=0.2) + geom_cladelab(node=filter(Annotation_phyla, Phylum_group=="Verrucomicrobia")$node, label= "Verrucomicrobia", angle=45, fontsize=3, vjust=0  ) + #Verrucomicrobia
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Lentisphaerae")$node  , fill=color_mapping["Lentisphaerae"], color=NA, alpha=0.2) + geom_cladelab(node=filter(Annotation_phyla, Phylum_group=="Lentisphaerae")$node, label= "Lentisphaerae", angle=55, fontsize=3, vjust=0  ) + #Lentisphaerae
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Candidatus_Melainabacteria")$node  , fill=color_mapping["Candidatus_Melainabacteria"], color=NA, alpha=0.2) + geom_cladelab(node=filter(Annotation_phyla, Phylum_group=="Candidatus_Melainabacteria")$node, label= "Candidatus_Melainabacteria", angle=0, fontsize=3, vjust=.1  ) +  #Candidatus_Melainabacteria
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Actinobacteria")$node  , fill=color_mapping["Actinobacteria"], color=NA, alpha=0.2) + geom_cladelab(node=filter(Annotation_phyla, Phylum_group=="Actinobacteria")$node, label= "Actinobacteria", angle=0, fontsize=3, vjust=.1  ) + #Actinobacteria
  new_scale_fill() + geom_tippoint(aes(col=Included_analysis ), size=0.75) + scale_color_manual( values = c( "FALSE"="grey", "TRUE" = "#E83845") ) + theme(legend.position="bottom")-> SGB_Tree2
  #geom_fruit( geom="geom_tile", mapping = aes(fill=Included_analysis), width=0.03,offset=1 ) + scale_fill_manual( values = c( "FALSE"="grey", "TRUE" = "#E83845") ) + 
  #new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Log10_number_samples), width=0.03,offset=0.05 ) + scale_fill_gradientn(colours = rev(c("darkred",  "orange", "white"))) #-> SGB_Tree



ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/SGB_tree.pdf",SGB_Tree)
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/SGB_tree_v2.pdf",SGB_Tree2)




read_tsv("/mnt/project/Error_iqtree", col_names=F) -> Check

tibble() -> Check_tibble
for (File in Check$X1){
  Path = paste0("/mnt/project/Phylogenies/", File, "/RAxML_bestTree.",File,".StrainPhlAn4.tre" )
  read.tree(Path) -> Tree
  Check_tibble %>% rbind( tibble(SGB=File, N = length(Tree$tip.label) ) ) -> Check_tibble
}
Check_tibble %>% arrange(desc(N))  





library(tidyverse)

load("/mnt/project/SGB_abundances/for_association.RData")

print(DF)
print(Metadata)

Data = left_join(DF, Metadata)
Result = tibble()
Results_prevalence = tibble()
for (Taxa in colnames(DF %>% select(-ID_anal) ) ){
  Data %>% select(Taxa, Continent) %>% group_by(Continent, !!sym(Taxa) ) %>% summarise(N = n())  -> D
  for (C in c("Africa", "Asia", "Europe", "North_America", "South_America", "Oceania") ){
    D %>% filter(Continent == C) -> Conti
    Conti[Conti[Taxa] ==1,]$N -> Y
    Conti[Conti[Taxa] ==0,]$N -> N
    Prevalence =  Y/(Y+N)
    tibble(Continent = C , Prevalence = Prevalence) %>% mutate(SGB = Taxa) %>% rbind(Results_prevalence, . ) -> Results_prevalence
    
  }
}  



#make histogram for panel B

Info %>% ggplot(aes(x= Number_samples_phylogeny)) + geom_histogram() + theme_bw() + scale_x_log10() + geom_vline(xintercept = 300, color= "red",linetype = "dashed", linewidth=2) + xlab("Number of samples per species-tree") -> histogram_all_samples
Info %>% filter(Number_samples_phylogeny >= 300) %>% ggplot(aes(x= Number_samples_phylogeny)) + geom_histogram() + theme_bw() + xlim(299, 23000) + scale_x_continuous(breaks = seq(300, max(Info$Number_samples_phylogeny), by = 1500))  + xlab("Number of samples per species-tree") + theme(axis.text.x = element_text(angle = 45, hjust = 1))  -> histogram_included_samples
histogram_all_samples / histogram_included_samples + plot_layout(widths = c(3, 1)) -> Histogram_plot
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/Figure2_histogramSamples.pdf", Histogram_plot)



####################################
##Heatmap associations#############
####################################


read_csv("/mnt/project/Make_Associations/Association/All_results.csv") -> Phenotype_assoc
read_csv("/mnt/project/Make_Associations/Association/All_results_filtered.csv") -> Phenotype_assoc_sign

read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypic_categories.tsv") -> Phenotype_categories

Phenotype_assoc %>% mutate(Ass_ID = paste0(SGB, "|", Phenotype )  ) -> Phenotype_assoc
Phenotype_assoc_sign %>% mutate(Ass_ID = paste0(SGB, "|", Phenotype )  ) -> Phenotype_assoc_sign

Phenotype_assoc %>% mutate(Phyl_evidence = ifelse(Ass_ID %in% Phenotype_assoc_sign$Ass_ID, T, F ) ) -> Phenotype_assoc
Phenotype_assoc %>% left_join(. , rename(Phenotype_categories, Phenotype=Phenos )) -> Phenotype_assoc

for (category in unique(Phenotype_assoc$Category) ){

Phenotype_assoc %>% filter(SGB %in% Phenotype_assoc_sign$SGB & Phenotype %in% Phenotype_assoc_sign$Phenotype ) %>% filter(Category == category) %>%
  ggplot(., aes(y = Phenotype, x = SGB, fill = Phyl_evidence))  + theme_bw() +
  geom_tile(color = "black") +
  #geom_text(aes(label = phylo_median), color = "white", size = 4) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=3), axis.text.y = element_text(size=3) ) + ggtitle(category) -> Plot
  
}


generate_plot <- function(category, Phenotype_assoc, Phenotype_assoc_sign) {
  plot <- Phenotype_assoc %>%
    filter(SGB %in% Phenotype_assoc_sign$SGB & Phenotype %in% Phenotype_assoc_sign$Phenotype) %>%
    filter(Category == category) %>%
    ggplot(aes(y = Phenotype, x = SGB, fill = Phyl_evidence)) +
    theme_bw() +
    geom_tile(color = "black", width = 1, height = 1) +
    coord_fixed() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
      axis.text.y = element_text(size = 5)
    ) + scale_fill_manual( values = c("grey", "#E83845") ) + 
    ggtitle(category)
  ggsave( paste0("/mnt/project/Make_Associations/Association/Results/Plots/Association_heatmap_",category, ".pdf"), plot )
  return(plot)
}
categories <- unique(Phenotype_assoc$Category)
plots <- lapply(categories, generate_plot, Phenotype_assoc = Phenotype_assoc, Phenotype_assoc_sign = Phenotype_assoc_sign)
combined_plot <- wrap_plots(plots, ncol = 1)  + plot_layout(guides = "collect")


##############################################  
##Get summary statistics per phenotype########
##############################################

Phenotype_assoc %>% filter(Covariates != "protein_intake") -> Phenotype_assoc  
unique(Phenotype_assoc$Phenotype)  %>% sort() -> Phenos
tibble(Phenos)


read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> AllPhenos
AllPhenos %>% select( c("study_name", "subject_id", "ID_anal", Phenos ) ) -> AllPhenos

#1. Filter out repeated
AllPhenos %>%  distinct(subject_id, .keep_all = TRUE)  %>%  distinct(ID_anal, .keep_all = TRUE) -> AllPhenos
#24,829 unique individuals


#2. Per phenotype: count how many studies is present and how many individiduals (prevalence)
Studies_with_pheno = tibble()
Prevalence_binary = tibble()
Stats_pheno = tibble()
for (Pheno in Phenos){
  AllPhenos %>% drop_na(Pheno) -> With_pheno
  With_pheno %>% group_by(study_name) %>% summarise(Number_individuals_with_pheno = n()) %>% mutate(Phenotype = Pheno)   -> N_per_study
  rbind(Studies_with_pheno, N_per_study) -> Studies_with_pheno
  With_pheno[Pheno] %>% as_vector() -> V
  if ( class(V) != "numeric" || length(unique(V)) < 4  ){
    #Get prevalence per study
    With_pheno %>% group_by(!!sym(Pheno), study_name ) %>% summarise(N = n()) %>% mutate(Phenotype = Pheno) %>% rename(Value = !!sym(Pheno)) %>%  ungroup()  -> Prevalence
    Prevalence <- Prevalence %>% left_join(N_per_study, by = c("study_name", "Phenotype" )) %>% mutate(Prevalence = N / Number_individuals_with_pheno)
    rbind(Prevalence_binary, Prevalence) -> Prevalence_binary
    
  } else {
    #Get summary stats per study
    With_pheno %>% group_by(study_name ) %>% summarise(Mean = mean(!!sym(Pheno)), SD = sd(!!sym(Pheno)), Max = max(!!sym(Pheno)), Min=min(!!sym(Pheno)) ) %>% mutate(Phenotype = Pheno) -> Stats
    Stats <- Stats %>% left_join(N_per_study, by = c("study_name", "Phenotype" )) 
    rbind(Stats_pheno, Stats) -> Stats_pheno
    
  }

}

Studies_with_pheno %>% write_tsv(. , "/mnt/project/Make_Associations/Association/Results/Info_Cohorts/Studies_with_phenotype.tsv")
Prevalence_binary %>% write_tsv(. , "/mnt/project/Make_Associations/Association/Results/Info_Cohorts/Prevalence_by_study.tsv")
Stats_pheno %>% write_tsv(. , "/mnt/project/Make_Associations/Association/Results/Info_Cohorts/SummaryStats_by_study.tsv")


############################
###Change names#############
############################

read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> AllPhenos
AllPhenos %>%  mutate( `ACE.Inhibitors` = coalesce( `ACE.Inhibitors`, `ACE.Inhibitor` )) -> AllPhenos
AllPhenos %>%  mutate( `Asthma` = coalesce( Asthma, asthma )) -> AllPhenos
AllPhenos %>%  mutate( `Hepatitis` = coalesce( Hepatitis, hepatitis )) -> AllPhenos
AllPhenos %>%  mutate( `Hypertension` = coalesce( Temp$Hypertension, Temp$hypertension )) -> AllPhenos
AllPhenos %>% rename(ankylosing_spondylitis = RA ) -> AllPhenos

AllPhenos %>% select(-c(ACE.Inhibitor, asthma, hepatitis, hypertension)) -> AllPhenos
AllPhenos %>% write_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv")

read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_tmp.tsv") -> Temp
Temp %>% filter(ID_anal %in% AllPhenos$ID_anal & subject_id %in%  AllPhenos$subject_id & id_used %in%  AllPhenos$id_used & subject_id %in%  AllPhenos$subject_id ) %>% distinct() -> Temp


coalesce( Temp$Hypertension, Temp$hypertension )





###################################
###Correlate geography matrices####
###################################
library(tidyverse)
library(ape)
library(vegan)
library(geosphere)

#Get coordinate info
#world <- map_data("world") %>% as_tibble() %>% filter(!region == "Antarctica")
cities = get('world.cities') %>% as_tibble() %>% filter(capital != 0) %>% rename(region=country.etc ) %>% rbind( tibble(name="Suva", region="Fiji", pop=NA,lat=-18.1405, long=178.44, capital=1 ))
read_csv("/mnt/project/Make_Associations/Phenotypes/Countries.csv", col_names = F) %>% rename(Country=X1, region=X2 ) -> Translation
Translation %>% rbind(tibble(Country = "KOR", region="Korea South" )) -> Translation
#Get tree
Clean_names =  function(Name){
  'Samples from 500FG FSK need everything after _ to be removed. If _metaphlan4 is attached at the end of the name, remove'
  if ( grepl("HV", Name)) {
    Name = str_split(Name, "_")[[1]][1]
  }   
  Name = str_replace(Name, "_metaphlan4", "") 
  return(Name)
}
read.tree("/mnt/project/Symlink_phylo/IQtree.t__SGB14546_group.TreeShrink.tre") -> Tree
if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }
Tree = midpoint.root(Tree)
Tree$tip.label = Tree$tip.label %>% sapply(Clean_names)
#Get data
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") %>% filter(! study_name == "ThomasAM_2019_c" ) -> Data
Data %>% filter(ID_anal %in% Tree$tip.label) -> Data
#remove repeated samples
Data %>% distinct(subject_id, .keep_all = T) %>% distinct(ID_anal, .keep_all = T)  -> Data2
#Merge data and geospatial information
Data2 %>% select(ID_anal, Country) %>% left_join( cities %>% left_join(Translation)  ) %>% group_by(ID_anal, region) %>% sample_n(1) -> Data3

#Merge data and tree information
Data3 %>% filter(ID_anal %in% Tree$tip.label) -> Data_anal
keep.tip(Tree, Data_anal$ID_anal) -> Tree
order_indices <- match(Tree$tip.label, as.character(Data_anal$ID_anal))
Data_anal[order_indices,] -> Data_anal

#calculate spatial and phylogenetic distances
Data_anal %>% ungroup() %>% select(long, lat) %>% as.matrix()  %>% distm(. ,  fun = distVincentySphere) -> distance_matrix
cophenetic.phylo(Tree) -> Distance_tree
#perform mantel test
mantel( as.dist(Distance_tree), as.dist(distance_matrix), method = "pearson", permutations=999, parallel=20 ) -> Test #0.3 r
Test$permutations
tibble(Rho= Test$statistic, P=Test$signif, Permutations= Test$permutations)




library(tidyverse)
library(ape)
library(vegan)
library(geosphere)

Run_analysis =  function(Tree_file, Threads, Perm =999 ){
  #Get coordinate info
  #world <- map_data("world") %>% as_tibble() %>% filter(!region == "Antarctica")
  cities = get('world.cities') %>% as_tibble() %>% filter(capital != 0) %>% rename(region=country.etc ) %>% rbind( tibble(name="Suva", region="Fiji", pop=NA,lat=-18.1405, long=178.44, capital=1 ))
  read_csv("/mnt/project//Make_Associations/Phenotypes/Countries.csv", col_names = F) %>% rename(Country=X1, region=X2 ) -> Translation
  Translation %>% rbind(tibble(Country = "KOR", region="Korea South" )) -> Translation
  #Get tree
  Clean_names =  function(Name){
    'Samples from 500FG FSK need everything after _ to be removed. If _metaphlan4 is attached at the end of the name, remove'
    if ( grepl("HV", Name)) { Name = str_split(Name, "_")[[1]][1] }   
    Name = str_replace(Name, "_metaphlan4", "") 
    return(Name)
  }
  read.tree(Tree_file) -> Tree
  if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }
  Tree = midpoint.root(Tree)
  Tree$tip.label = Tree$tip.label %>% sapply(Clean_names)
  #Get data
  read_tsv("/mnt/project//Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") %>% filter(! study_name == "ThomasAM_2019_c" ) -> Data
  Data %>% filter(ID_anal %in% Tree$tip.label) -> Data
  #remove repeated samples
  Data %>% distinct(subject_id, .keep_all = T) %>% distinct(ID_anal, .keep_all = T)  -> Data2
  #Merge data and geospatial information
  Data2 %>% select(ID_anal, Country, Continent, study_name) %>% left_join( cities %>% left_join(Translation)  ) %>% group_by(ID_anal, region) %>% sample_n(1) -> Data3
  
  #Merge data and tree information
  Data3 %>% filter(ID_anal %in% Tree$tip.label) -> Data_anal
  keep.tip(Tree, Data_anal$ID_anal) -> Tree
  order_indices <- match(Tree$tip.label, as.character(Data_anal$ID_anal))
  Data_anal[order_indices,] -> Data_anal
  
  ggtree(Tree, layout = "circular")  %<+% Data_anal + geom_tippoint(aes(col=Country), size=0.5 ) + scale_color_manual( values = C_colors  ) + 
    new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.01,offset=0.1 ) + scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30"))
  
  #calculate spatial and phylogenetic distances
  #Data_anal %>% ungroup() %>% select(long, lat) %>% as.matrix()  %>% distm(. ,  fun = distVincentySphere) -> distance_matrix
  #cophenetic.phylo(Tree) -> Distance_tree
  #perform mantel test
  #mantel( as.dist(Distance_tree), as.dist(distance_matrix), method = "pearson", permutations=Perm, parallel=Threads ) -> Test #0.3 r
  #tibble(Rho= Test$statistic, P=Test$signif, Permutations= Test$permutations) -> Results
  return( Country_plot )
  
  
}

read_tsv("/mnt/project/Make_Associations/Association/Tree_sizes2.txt") %>% filter(! N_tips < 300  ) -> Trees_check


SGB_from_tree = function(filename){
  split_filename <- unlist(strsplit(filename, "\\."))
  # Extract the second element
  extracted_string <- split_filename[2]
}
Files <- list.files(path = "/mnt/project/Symlink_phylo/", pattern = "\\.tre$", full.names = TRUE)
Results = tibble()
for (Tree_file in Files){
  if (! basename(Tree_file) %in% Trees_check$Tree ){ next }
  SGB = SGB_from_tree(Tree_file)
  print(SGB)
  Out = paste0("/mnt/project/Make_Associations/Association/Results/Geography_cor/", SGB, ".pdf")
  #if(file.exists(Out)){ print("Already done") ; next }
  Run_analysis(Tree_file, Threads = 20, Perm = 2000) -> Result #%>% mutate(SGB=SGB) -> Result
  ggsave(Out, Result)
  #Result %>% rbind(Results, . ) -> Results
}
write_tsv(Results, "/mnt/project/Make_Associations/Association/Results/Geography_cor/All_cor.tsv")

Files <- list.files(path = "/mnt/project/Make_Associations/Association/Results/Geography_cor/", pattern = "\\.tsv$", full.names = TRUE)
Results_merge = tibble()
for (File in Files){
  if (grepl("All_cor", File)){ next }
  if (grepl("merged", File)){ next }
  read_tsv(File, show_col_types = FALSE) -> Fi
  if (dim(Fi)[2] < 5){ print(File) ; next } 
  Fi %>% rbind(Results_merge, .) -> Results_merge
}
Results_merge %>% ggplot(aes(x=Rho, y=Rho_spear)) + geom_point() + geom_abline() + theme_bw()
Results_merge %>% mutate(FDR = p.adjust(P, "fdr"), FDR_spear = p.adjust(P_spear, "fdr") ) -> Results_merge
Results_merge %>% write_tsv("/mnt/project/Make_Associations/Association/Results/Geography_cor/All_cor_withspear.tsv")


c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
read_tsv("/mnt/project//Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") %>% filter(! study_name == "ThomasAM_2019_c" ) -> Data
unique(Data$Country) -> C_name
c25[1:length(C_name)] -> C_colors
names(C_colors) = C_name


######################################################
##########Get normalized pairwise distances###########
######################################################
library(ape)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 ){ print("Please pass a SGB and tree file in the command line: Rscript scripts/Extract_norm_distance.R {SGB} {Tree}") ; q() }
Tree = read.tree( args[1] )
SGB = args[2]
if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }

Items = Tree$tip.label
Distance_matrix = cophenetic.phylo(Tree)
upper_triangular <- Distance_matrix[upper.tri(Distance_matrix, diag = TRUE)]
# Get row and column indices for upper triangular elements
row_indices <- row(Distance_matrix)[upper.tri(Distance_matrix, diag = TRUE)]
col_indices <- col(Distance_matrix)[upper.tri(Distance_matrix, diag = TRUE)]
# Create a data frame with row indices, column indices, and distances
pairwise_distances <- tibble( Row = rownames(Distance_matrix)[row_indices], Column = colnames(Distance_matrix)[col_indices], Distance = upper_triangular)
# Print the resulting data frame
pairwise_distances %>% filter(! Row == Column) -> pairwise_distances
Median_pairwise = median(pairwise_distances$Distance)
pairwise_distances %>% mutate(Distance = Distance/Median_pairwise ) -> Normalized_distance
Normalized_distance %>% write_tsv("Normalized_distances/",SGB,".tsv", col_names = F)


#######################################################
#################Format realtionship data#############
#######################################################

read_tsv("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Strain_sharedness/Data/Info_relation_raw.tsv")
Info = read_tsv("/mnt/project/Strain_sharedness/Data/Info_relation_raw.tsv")
Info_dmp = read_tsv("/mnt/project/Strain_sharedness/Data/cohousing_dmp2.tsv")
Info_vallescolomer = read_tsv("/mnt/project/Strain_sharedness/Data/Metadata_VallesColomer.tsv")
Info_vallescolomer %>% rename(ID_anal = sampleID, family=familyID) %>% select(ID_anal, subjectID, family) -> Info_vallescolomer2

Info  %>% drop_na() %>% group_by(family) %>% summarise(n())
Info <- Info %>% left_join(Info_dmp, by = "ID_anal") %>%
  mutate(family = coalesce(family.y, family.x)) %>%
  select(colnames(Info))

Info %>% left_join(Info_vallescolomer2 %>% select(-subjectID), by = "ID_anal") %>%
  mutate(family = coalesce(family.y, family.x)) %>%
  select(colnames(Info)) -> Info

#Info %>% write_tsv("/mnt/project/Strain_sharedness/Data/Information_complete.tsv")
#Info %>% select(ID_anal,subject_id, family, study_name) %>% write_tsv("/mnt/project/Strain_sharedness/Data/Input_intrapopulation.tsv")


read_tsv("/mnt/project/Strain_sharedness/Data/mdata_transm_27datasets_stool_healthy_forscript_HOUSEHOLD_3M_20200729.tsv") -> data_mireia

data_mireia %>% rename(ID_anal =  sampleID) %>% select(ID_anal, role, time_point) %>% left_join(Info, . ) %>% mutate( role = ifelse(is.na(role),"U",role) ) -> Info
write_tsv(Info, "/mnt/project/Strain_sharedness/Data/Information_complete.tsv")

Info %>% mutate(family.1=family) %>% select(ID_anal, study_name,subject_id, family, family.1, role, time_point) %>% write_tsv("/mnt/project/Strain_sharedness/Data/Input_intrapopulation.tsv")


##############################################################
#################Country/Contient median distance#############
##############################################################

library(data.table)
# Read Distance as data.table
Distance_file = "/mnt/project/Strain_sharedness/Normalized_distances/t__SGB10068.tsv"
SGB <- basename(Distance_file) %>% tools::file_path_sans_ext(.)
Distance <- fread(Distance_file)
# Read AllPhenos as data.table
AllPhenos <- fread("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv")

# Create a set of unique IDs from Distance$X1 and Distance$X2
unique_ids <- unique(c(Distance$V1, Distance$V2))

# Filter Phenos_unique using unique_ids and convert to data.table
Phenos_unique = AllPhenos %>% filter( ID_anal %in% unique_ids ) %>% distinct(sample_id, .keep_all = T) %>% distinct(ID_anal, .keep_all = T) %>% select(c(ID_anal,Country, Continent) )
Phenos_unique = Phenos_unique %>% data.table() %>% as.tibble() %>% as.data.table() 
# Set key for fast joining
setkey(Distance, V1)
setkey(Phenos_unique, ID_anal)

# Perform left join with Distance and Phenos_unique for X1
Distance <- Distance[Phenos_unique, nomatch = 0]

# Set key for fast joining
setkey(Distance, V2)
setkey(Phenos_unique, ID_anal)

# Perform left join with Distance and Phenos_unique for X2
Distance <- Distance[Phenos_unique, nomatch = 0]

# Rename columns and select required columns
Distance <- Distance[, .(V1, V2, V3 ,Country_v1=Country, Continent_v1=Continent, Country_v2=i.Country, Continent_v2=i.Continent)]

#1. Get mean/median per country/Continent
Distance[, Continent_ID := paste(pmin(Continent_v1, Continent_v2), pmax(Contienent_v1, Continent_v2), sep = "-")]
Distance[, Country_ID := paste(pmin(Country_v1, Country_v2), pmax(Country_v1, Country_v2), sep = "-")]

Distance %>% group_by(Contient_ID) %>% sumamrise(median(V3))

result_continents <- Distance[, .(Median_distance_by_continent = median(V3)), by = Continent_ID]
result_countries <- Distance[, .(Median_distance_by_country = median(V3)), by = Country_ID]


Output = paste0("/mnt/project/Strain_sharedness/Distance_countryNcontinent/",SGB, "_continent.tsv")
write_tsv(result_continents, Output)

Output = paste0("/mnt/project/Strain_sharedness/Distance_countryNcontinent/",SGB, "_country.tsv")
write_tsv(result_countries, Output)




#################################################################################################
#1. Get a dendogram from median phylogenetic distances by performing a hierarchical clustering###
#################################################################################################
Get_dendogram = function(result_df){
  #Make sure always same names
  colnames(result_df) = c("Continent_ID", "Median_distance_by_continent")
  #Generate two columns
  split_df <- separate(result_df, Continent_ID, into = c("Continent1", "Continent2"), sep = "-")
  #Make a wide format
  wide_df =pivot_wider(split_df, names_from = "Continent2", values_from = "Median_distance_by_continent")
  wide_df = as.data.table(wide_df)
  #Make matrix
  distance_matrix = as.matrix(wide_df, rownames="Continent1")
  row_indices <- match( sort(rownames(distance_matrix)), rownames(distance_matrix))
  col_indices <- match( sort(rownames(distance_matrix)), colnames(distance_matrix))
  distance_matrix <- distance_matrix[row_indices, col_indices]
  #Make lowertrinagle to upper triangle
  distance_matrix[lower.tri(distance_matrix)] <- t(distance_matrix)[lower.tri(distance_matrix)]
  #Plot and generate dendogram
  pheatmap::pheatmap(distance_matrix)
  dendogram = hclust(as.dist(distance_matrix))
  plot(dendogram)
  dendogram = as.phylo(dendogram)
  return(dendogram)
}


###########################################################################################################
#2. Compare distances to different continents in western/non-western subpopulations #######################
###########################################################################################################
library(gghalves)
library(data.table)

# Read Distance as data.table
Distance_file = "/mnt/project/Strain_sharedness/Normalized_distances/t__SGB10068.tsv"
SGB <- basename(Distance_file) %>% tools::file_path_sans_ext(.)
Distance <- fread(Distance_file)
# Read AllPhenos as data.table
AllPhenos <- fread("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv")
#Add westernization info
VallesColomer = fread("/mnt/project/Strain_sharedness/Data/Westernization/VallesColomerM_2023.tsv") %>% filter(subcohort == "VallesColomerM_2023_GNB")
KurK = fread("/mnt/project/Strain_sharedness/Data/Westernization/KurK_2020.tsv")
LiuW = fread("/mnt/project/Strain_sharedness/Data/Westernization/LiuW_2016.tsv")

rbind( (VallesColomer %>% rename(ID_anal = sample_id) %>% select(ID_anal, non_westernized) ),   (KurK %>% rename(ID_anal = sample_id) %>% select(ID_anal, non_westernized))  )  %>%
  rbind(  (LiuW %>% rename(ID_anal = sample_id) %>% select(ID_anal, non_westernized))) -> NonWestern

AllPhenos  %>% left_join( NonWestern, by="ID_anal") %>% group_by(study_name, non_westernized) -> AllPhenos


# Create a set of unique IDs from Distance$X1 and Distance$X2
unique_ids <- unique(c(Distance$V1, Distance$V2))

# Filter Phenos_unique using unique_ids and convert to data.table
Phenos_unique = AllPhenos %>% filter( ID_anal %in% unique_ids ) %>% distinct(sample_id, .keep_all = T) %>% distinct(ID_anal, .keep_all = T) %>% select(c(ID_anal,Country, Continent, non_westernized) )
Phenos_unique = Phenos_unique %>% data.table() %>% as.tibble() %>% as.data.table() 
# Set key for fast joining
setkey(Distance, V1)
setkey(Phenos_unique, ID_anal)

# Perform left join with Distance and Phenos_unique for X1
Distance <- Distance[Phenos_unique, nomatch = 0]

# Set key for fast joining
setkey(Distance, V2)
setkey(Phenos_unique, ID_anal)

# Perform left join with Distance and Phenos_unique for X2
Distance <- Distance[Phenos_unique, nomatch = 0]

# Rename columns and select required columns
Distance <- Distance[, .(V1, V2, V3 ,Country_v1=Country, Continent_v1=Continent,non_westernized_v1=non_westernized ,Country_v2=i.Country, Continent_v2=i.Continent,non_westernized_v2=i.non_westernized) ]


#Per each of the datasets. If we compare western to non-western, do we see non-western closer to distances of western?
##LIUW##
Distance %>% as_tibble() %>% dplyr::filter(V1 %in%  LiuW$sample_id |  V2 %in% LiuW$sample_id  ) -> LiuW_distances

LiuW_distances$non_westernized_v1[is.na(LiuW_distances$non_westernized_v1)] <- LiuW_distances$Continent_v1[is.na(LiuW_distances$non_westernized_v1)]
LiuW_distances$non_westernized_v2[is.na(LiuW_distances$non_westernized_v2)] <- LiuW_distances$Continent_v2[is.na(LiuW_distances$non_westernized_v2)]

LiuW_distances %>% mutate(Comparison = paste(pmin(non_westernized_v1, non_westernized_v2), pmax(non_westernized_v1, non_westernized_v2), sep = "-")) -> LiuW_distances
LiuW_distances %>% filter(grepl("North_America|Europe|Asia", Comparison) | Comparison %in% c("no-yes")  ) %>% ggplot(aes(y=V3, x=Comparison )) + 
  geom_half_boxplot() + theme_bw() + coord_flip() +
  gghalves::geom_half_violin(side = "t", aes(fill = Comparison), width = 0.5) 
#LiuW_distances %>% filter(grepl("Nort_America|Europe|Asia", Comparison) | Comparison %in% c("no-yes")  )  %>%
#  lm( V3 ~  )

##KurK##
Distance %>% as_tibble() %>% dplyr::filter(V1 %in%  KurK$sample_id |  V2 %in% KurK$sample_id | Country_v1=="IND" |  Country_v2=="IND" ) -> KurK_distances
KurK_distances %>% mutate(non_westernized_v1 = ifelse(is.na(non_westernized_v1) & Country_v1=="IND", Country_v1, non_westernized_v1) ) %>%
  mutate(non_westernized_v2 = ifelse(is.na(non_westernized_v2) & Country_v2=="IND", Country_v2, non_westernized_v2)) -> KurK_distances
KurK_distances$non_westernized_v1[is.na(KurK_distances$non_westernized_v1)] <- KurK_distances$Continent_v1[is.na(KurK_distances$non_westernized_v1)]
KurK_distances$non_westernized_v2[is.na(KurK_distances$non_westernized_v2)] <- KurK_distances$Continent_v2[is.na(KurK_distances$non_westernized_v2)]

KurK_distances %>% mutate(Comparison = paste(pmin(non_westernized_v1, non_westernized_v2), pmax(non_westernized_v1, non_westernized_v2), sep = "-")) -> KurK_distances
KurK_distances %>% filter(grepl("North_America|Europe|Asia|India", Comparison) ) %>% ggplot(aes(y=V3, x=Comparison )) + 
  geom_half_boxplot() + theme_bw() + coord_flip() +
  gghalves::geom_half_violin(side = "t", aes(fill = Comparison), width = 0.5) 


###############################################
#####Check longitudinal samples in missing#####
###############################################
Missing = read_tsv("/mnt/project/Strain_sharedness/Missing_thresholds.txt", col_names = F)

Check_longitudinal = function(SGB){
  print(SGB)
  matching_files <- list.files("/mnt/project/Symlink_phylo/", pattern = paste0(".",SGB,".Tree"), full.names = TRUE)
  read.tree(matching_files) -> Tree
  if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }
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
  
  read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> AllPhenos
  AllPhenos %>% filter(ID_anal %in% Tree$tip.label ) -> Phenos
  
  result <- tryCatch(
    {
      Phenos %>%
        group_by(sample_id) %>%
        summarise(N = n()) %>%
        print() %>%
        filter(!is.na(sample_id)) %>%
        filter(N > 1) %>%
        mutate(SGB = SGB)
    },
    error = function(e) {
      tibble(SGB = SGB, N = NA)
    }
    
  )
  Phenos %>% filter(ID_anal %in% result$sample_id) %>% group_by(study_name) %>% summarise(n())
  
  
  return(result)
}

library(doParallel)
library(foreach)
library(ape)
cl <- makeCluster(10)
registerDoParallel(cl)

# Create a foreach loop with %dopar%
result_list <- foreach(SGB = Missing$X1, .combine = c, .packages = c("tidyverse", "ape")  ) %dopar% {
  cat(SGB)
  Check_longitudinal(SGB)
}

# Stop the parallel backend
stopCluster(cl)
registerDoSEQ()  # Return to sequential processing

# Combine the results from different cores
combined_result <- unlist(result_list)

Missing$X1 %>% sapply(Check_longitudinal) -> Result


read_tsv("/mnt/project/Strain_sharedness/Longitudinal_samples_missing.tsv") -> D
D %>% group_by(N, SGB) %>% summarise(N2 = n()) %>% arrange(desc(N2)) %>% filter(N2>50)
#For each of them, check if timepoints are within 6 months
# t__SGB1814 two datasets LeChatelierE_2013/NielsenHB_2014 with same participants
#t__SGB1815 same
#t__SGB1934 same
#t__SGB1871 same




####For each missing SGB, extract distribution of distances for non-longitudinal, take 3rd percentile


library(ape)
library(tidyverse)

SGB = "t__SGB1871"
Tree = "/mnt/project/Symlink_phylo//RAxML_bestTree.t__SGB1871.TreeShrink.tre"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 ){ print("Please pass a SGB and tree file in the command line: Rscript scripts/Extract_norm_distance.R {SGB} {Tree}") ; q() }
Tree = read.tree( args[1] )
SGB = args[2]
if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }
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
#Read meta
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> AllPhenos
AllPhenos %>% filter(ID_anal %in% Tree$tip.label ) -> Phenos
Phenos %>% distinct(ID_anal, .keep_all = T) %>% distinct(sample_id, .keep_all = T) -> Phenos
keep.tip(Tree, Phenos$ID_anal) -> Tree

Items = Tree$tip.label
Distance_matrix = cophenetic.phylo(Tree)
Distance_matrix[upper.tri(Distance_matrix, diag = F)] %>% as.vector() -> Distance
Distance 
ggplot() + geom_histogram(aes(x=Distance)) + theme_bw()
third_percentile <- quantile(Distance, probs = 0.03)
return(tibble(SGB=SGB, Threshold = third_percentile))




read.table("/mnt/project/Strain_sharedness/Human_phylogeny/Data/Kinship.tsv", header = T,sep = "\t") -> Distance
pheatmap::pheatmap(Distance)


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

Plot_enrichment = plotEnrichment(taxa_list[["f__Lachnospiraceae"]], ranks) + labs(title="Lachnospiraceae") + theme( text = element_text(size = 21),  plot.title = element_text(size = 26), axis.title.x = element_text(size = 21),  axis.title.y = element_text(size = 21)  ) 
plotEnrichment(taxa_list[["o__Clostridiales"]], ranks) + labs(title="Clostridiales")
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/Lachnospiraceae_enrichment_geography.pdf",Plot_enrichment)


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

##enrichments in subanalyses?
Overall_associations =  read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Europe.csv" ) %>% mutate(ID = paste0(Phenotype,SGB) )
Significant_associations = read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Europe_filtered.csv" )%>% mutate(ID = paste0(Phenotype,SGB) )
Overall_associations %>% mutate(FDR = ifelse(ID  %in% Significant_associations$ID, 0, 1 )) -> Overall_associations






#######################################
#### Proportion transmission events ###
#######################################

file_list <- list.files(path = "/mnt/project/Strain_sharedness/Transmission", pattern = "Annotation_info\\.tsv$", recursive = TRUE, full.names = TRUE)

File = file_list[1]
SGB = basename(dirname(File))

read_tsv(File) -> Info
Info %>% group_by(are_mother_infant) %>% summarise(N = n()) %>% filter(are_mother_infant == T  )
Info %>% group_by(are_families) %>% summarise(N = n()) %>% filter(are_families == T  )
Info %>% group_by(are_intrastudy) %>% summarise(N = n()) %>% filter(are_intrastudy == T  )
Info %>% group_by(are_interstudy) %>% summarise(N = n()) %>% filter(are_interstudy == T  )


######################################
#### Sumamry assocations per group ###
######################################

read_csv("/mnt/project/Make_Associations/Association/All_results_filtered.csv") -> Sign
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypic_categories.tsv") -> Cat

Sign %>%
  left_join(Cat %>% rename(Phenotype = Phenos)) %>%
  group_by(Phenotype) %>%
  summarise(N_associations = n(), Category) %>%
  distinct(.keep_all = TRUE) %>%
  arrange(desc(N_associations)) %>%
  head(10) %>%
  ggplot(aes(y = N_associations, x = reorder(Phenotype, -N_associations), fill = Category)) +
  geom_bar(stat = "identity") +
  theme_bw() + coord_flip() + scale_fill_manual(values = c("#80678e", "#db8c6d", "#ef7b8e" ) ) + ylab("Number of associated phylogenies") + xlab("Human Phenotype")


###Format input strain sharedness#########

read_tsv("/mnt/project/Strain_sharedness/Data/Information_complete.tsv") -> Info 
Info %>% mutate(family = ifelse(is.na(family), ID_anal, family ) ) -> Info
Info %>% group_by(subject_id) %>% mutate(time_point = row_number()) %>% ungroup() -> Info
write_tsv(Info, "/mnt/project/Strain_sharedness/Data/Information_complete.tsv")



###############################
####Follow=up diseaseses######
##############################

Disease = "IBD"
Metadata = read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv")
Metadata %>% filter( !!sym(Disease) == 1 ) %>% group_by(study_name) %>% summarise(N = n()) -> studies
Metadata %>% filter(study_name %in% studies$study_name) %>% filter(!is.na( !!sym(Disease) )) %>%  group_by(study_name) %>% summarize(proportion = mean(!!sym(Disease))) -> proportions

Studies_keep = filter( proportions, proportion>0)$study_name

Metadata %>% filter(study_name %in% Studies_keep) %>% filter(!is.na( !!sym(Disease) )) -> Metadata

Metadata %>% group_by( study_name,Country ) %>% summarise(n()) -> Countries

Metadata %>% distinct(study_name, .keep_all = T) %>% group_by(Continent) %>% summarise(n()) -> Continents
Metadata %>% distinct(study_name, .keep_all = T) %>% group_by(Country) %>% summarise(n(), Continent) %>% distinct(.keep_all=T) -> Countries
Continents 


###Plot geography assocaitions####
#Show distribution and significance
#Add distribution of sharedness somehow
#show some tree examples

library(tidyverse)
library(ape)
library(phytools)

set.seed(914513)

Remove_repeated = function(Metadata, pruned_tree){
  'Check which samples are repeated (either duplicates or longitudinally) and only picks one from the repeated'
  set.seed(89777)
  Metadata %>% group_by(subject_id) %>% summarise(N = n()) %>% arrange(desc(N)) %>% filter(N > 1 ) -> Replicates
  remove = tibble()
  for (i in Replicates$subject_id){
    Metadata %>% filter(subject_id == i) -> Filtered ; Samples = Filtered$ID_anal
    Keep = sample(Samples, 1)
    remove = c(remove, Samples[!Samples==Keep] ) 
  }
  remove = unlist(remove)
  Metadata %>% filter(! ID_anal %in% remove) -> Metadata
  pruned_tree %>% drop.tip(remove) -> pruned_tree
  return(list(Metadata, pruned_tree))
}


Geography_summary = "/mnt/project/Make_Associations/Association/Results/Geography_cor/Association_results_merged.tsv"
read_tsv(Geography_summary) -> Geography_summary
Info_taxa =  "/mnt/project/Make_Associations/Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2"
Trees = read_tsv("/mnt/project/Make_Associations/Association/Tree_sizes2.txt") %>% 
  mutate(SGB = str_split(Tree, "\\.", simplify = TRUE)[, 2]) #%>% mutate(SGB= gsub("t__", "", SGB))
Metadata = read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv")
Metadata %>% filter(! study_name %in% c("ThomasAM_2019_c","500FG_FSK")  ) %>% group_by(ID_anal) %>% sample_n(1) %>% distinct(ID_anal, .keep_all = T) %>% ungroup()  -> Metadata

Proportion_sharedness =  read_tsv("/mnt/project/Strain_sharedness/Tranmission_rates.tsv")


Geography_summary %>% left_join(Trees) %>% left_join(Proportion_sharedness %>% filter(Intracountry_diffStudy==T ) ) -> Geography_summary
#histogram distribution of correlations
Geography_summary %>% ggplot() + geom_histogram(aes(x=Rho, fill=FDR<0.05 ), alpha=0.5, position = "identity" ) + theme_bw()  + scale_fill_manual( values = c("grey", "#E83845") ) + 
  geom_bar(aes(x=Rho, y=Shared) ,stat = "summary_bin", fun = mean)

Info$N_tips

Geography_summary %>% mutate(Rho_Bin = cut(Rho, breaks = 21, labels = FALSE)) %>% group_by(Rho_Bin, FDR<0.05) %>% summarise(Frequency = n(),Shared_bin = mean(Shared), Rho) -> Input_plot
C = 0.001

Input_plot %>% group_by(Rho_Bin, `FDR < 0.05`) %>%
  summarize(
    Frequency = first(Frequency),
    Shared_bin = first(Shared_bin),
    Rho = mean(Rho)
  )  %>% ggplot(aes(x=Rho)) + geom_bar(aes(y=Frequency, fill = `FDR < 0.05` ), color="black", stat="identity", position="identity", width=0.035, alpha=0.7) +
  geom_point(aes(y=Shared_bin/C,shape=`FDR < 0.05`),stat = "identity") + geom_line(aes(y=Shared_bin/C,shape=`FDR < 0.05`),stat = "identity") +
  scale_y_continuous( name = "Frequency Rho", sec.axis = sec_axis(~.*C, name="Mean sharedness")) +
  theme_bw() + scale_fill_manual( values = c("grey", "#E83845") ) +
  theme(
    text = element_text(size = 16),  # Set the text size for all elements
    plot.title = element_text(size = 20, face = "bold"),  # Set title text size
    axis.title = element_text(size = 18, face = "bold"),  # Set axis title text size
    axis.text = element_text(size = 14),  # Set axis text size
    axis.ticks = element_line(size = 1.5),  # Set axis tick size
    legend.text = element_text(size = 14),  # Set legend text size
    legend.title = element_text(size = 16, face = "bold")  # Set legend title text size
  ) -> Histogram



#Show example plots
read_tsv(Info_taxa, col_names =F) %>% rename(SGB=X1, Taxonomy=X2) -> Info_taxa
Geography_summary %>% mutate(SGB = gsub("t__", "", SGB)) %>% left_join(Info_taxa) %>% left_join(Trees) -> Geography_summary
Geography_summary %>% arrange(desc(Rho)) %>% select(Rho, SGB, Taxonomy)

To_plot= c("SGB1408", "SGB15145")
S = To_plot[1]

Make_tree = function(Geography_summary, S, Metadata, tip=1.4){
  Geography_summary %>% filter(SGB %in%  c(S, paste0("t__",S) ) ) -> Filtered
  Tree = read.tree(paste0("/mnt/project/Symlink_phylo/", Filtered$Tree))
  if (class(Tree) == "multiPhylo" ){ Tree = Tree[[1]] }
  #1. Clean names
  Tree = midpoint.root(Tree)
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
  #2. Remove repeated
  Metadata %>% filter(ID_anal %in% Tree$tip.label ) %>% filter(!Country == "NLD") -> Metadata_anal
  pruned_tree =  keep.tip(Tree, Metadata_anal$ID_anal)
  Unrepeated = Remove_repeated(Metadata_anal, pruned_tree )
  Metadata_anal = Unrepeated[[1]] ; pruned_tree = Unrepeated[[2]]
  Metadata_anal %>% select(ID_anal, Country, Continent) %>% rename(sample_id = ID_anal) %>% drop_na() -> Metadata_anal
  pruned_tree %>% keep.tip(Metadata_anal$sample_id ) -> pruned_tree
  #3. Plot
  ggtree(pruned_tree, layout="fan", open.angle=15, size=0.1) %<+% Metadata_anal -> p
  p +  geom_tiplab(aes(col=Continent, label=Country ), size=1.4) + scale_color_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30") )
}
     
SGB1408_tree = Make_tree(Geography_summary, "SGB1408", Metadata ) #rho: 0.641
SGB15145_tree = Make_tree(Geography_summary, "SGB15145", Metadata ) #rho: 0.609
t__SGB4552_group_tree = Make_tree(Geography_summary, "t__SGB4552_group", Metadata )  #1 t__SGB4552_group 
t__SGB4041_tree = Make_tree(Geography_summary, "t__SGB4041", Metadata )  #2 t__SGB4041 

Torques1_tree = Make_tree(Geography_summary, "SGB4563_group", Metadata, tip=1 ) #rho: 0.379
Rectales1_tree = Make_tree(Geography_summary, "SGB4933_group", Metadata, tip=1 ) #rho: 0.3
Aerofaciens1_tree = Make_tree(Geography_summary, "SGB14546_group", Metadata, tip=1 ) #aerofaciens, highligted in Suzuki Science paper

Ventrisium = Make_tree(Geography_summary, "SGB5045", Metadata, tip=1 ) #Top one with > 3K samples
Lachnospiraceae_unc = Make_tree(Geography_summary, "SGB4910", Metadata, tip=1 ) #Top one with > 6K samples

Onderdonkii = Make_tree(Geography_summary, "SGB2303", Metadata, tip=1 )  #SGB2303/Alistipes_onderdonkii, low geogrpahical effect (0.0242), sign FDR 0.000790, over 10K samples (12042)

#In ordre: Troques (clear caldes),  rectale (previously seen in Nicola), aerofaciens (highlighted in Suzuki), Ventrisium (large effect in more than 3K), Lachnospiraceae (top oen with >6K), Alistipes_onderdonkii (low one with over 10K samples)
SGB_geography_highlight  = c("SGB4563_group", "SGB4933_group", "SGB14546_group", "SGB5045", "SGB4910", "SGB2303" ) 


#Save plot
Prefix = "/mnt/project/Make_Associations/Association/Results/General_plots/Figure_geography_effect/"
ggsave( paste0(Prefix, "Histogram_rho.pdf"), Histogram )

ggsave( paste0(Prefix, "Tree_SGB1408.pdf"), SGB1408_tree )
ggsave( paste0(Prefix, "Tree_SGB15145.pdf"), SGB15145_tree )
ggsave(paste0(Prefix, "Tree_SGB4563_group.pdf"), Torques1_tree )

ggsave(paste0(Prefix, "Tree_SGB4933_group.pdf"), Rectales1_tree )
ggsave(paste0(Prefix, "Tree_SGB14546_group.pdf"), Aerofaciens1_tree ) 


ggsave(paste0(Prefix, "Tree_SGB5045.pdf"), Ventrisium )
ggsave(paste0(Prefix, "Tree_SGB4910.pdf"), Lachnospiraceae_unc ) 
ggsave(paste0(Prefix, "Tree_SGB2303.pdf"), Onderdonkii ) 



###Summaries associations
Significant_associations = read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/All_results_filtered.csv" )%>% mutate(ID = paste0(Phenotype,SGB) )
Pheno_info = read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypic_categories.tsv")
Significant_associations %>% left_join( Pheno_info%>% rename(Phenotype=Phenos) ) -> Significant_associations
Significant_associations%>% group_by(Category) %>% summarise(n())
Significant_associations %>% mutate(ID = paste0(SGB, "_", Phenotype) ) -> Significant_associations


Significant_associations %>% filter(Category == "Biochemical") %>% group_by(Phenotype) %>% summarise(n())
Significant_associations %>% filter(Category == "Lifestyle_and_exposome") %>% group_by(Phenotype) %>% summarise(n())

SGB_taxonomy %>% mutate(SGB = paste0("t__", SGB)) %>% filter(SGB %in% filter(Significant_associations, Phenotype %in% c("NO2") )$SGB  ) %>% select(Taxonomy)
SGB_taxonomy %>% mutate(SGB = paste0("t__", SGB)) %>% filter(SGB %in% c("t__SGB15271",  "t__SGB14042", "t__SGB15090" ) ) %>% select(SGB, Taxonomy)


Significant_associations %>% filter(Phenotype %in% c("Calprotectin", "triglycerides", "bilirubin", "albumine") )


Significant_associations %>% filter(Category == "Other") %>% group_by(Phenotype) %>% summarise(n())
Significant_associations %>% filter(Phenotype %in% c("BristolFreq", "BristolType") ) %>% select(ID)


Significant_associations %>% filter(Phenotype == "BirthMode")


###Check the disease-model
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_filtered.csv") -> Disease_model
Disease_model %>% group_by(Phenotype) %>% summarise(n())
Disease_model %>% mutate(SGB = str_replace(SGB, "t__", "") ) %>% left_join(SGB_taxonomy ) -> Disease_model


Plot_heatmap = function(Disease_model){
Disease_model %>%
    ggplot(aes(y = Phenotype, x = species, fill = elpd_diff)) +
    theme_bw() +
    geom_tile(color = "black", width = 1, height = 1) +
    coord_fixed() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5)
    ) 

}


read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_Europe_filtered.csv") %>%  mutate(SGB= str_replace(SGB, "t__", "") ) %>% left_join(SGB_taxonomy %>% select(SGB, species, Taxonomy)  )
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_NorthAmerica_filtered.csv") %>%  mutate(SGB= str_replace(SGB, "t__", "") ) %>% left_join(SGB_taxonomy %>% select(SGB, species, Taxonomy)  )
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_Asia_filtered.csv") %>%  mutate(SGB= str_replace(SGB, "t__", "") ) %>% left_join(SGB_taxonomy %>% select(SGB, species, Taxonomy)  )

read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disese_perContinent/Results_PerContinent_Europe_filtered.csv") %>%  mutate(SGB= str_replace(SGB, "t__", "") ) %>% left_join(SGB_taxonomy %>% select(SGB, species, Taxonomy)  ) -> EuropeAss
Plot_heatmap(EuropeAss)
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disese_perContinent/Results_PerContinent_Asia_filtered.csv") %>%  mutate(SGB= str_replace(SGB, "t__", "") ) %>% left_join(SGB_taxonomy %>% select(SGB, species, Taxonomy)  ) -> AsiaAss
Plot_heatmap(AsiaAss)
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disese_perContinent/Results_PerContinent_NorthAmerica_filtered.csv") %>%  mutate(SGB= str_replace(SGB, "t__", "") ) %>% left_join(SGB_taxonomy %>% select(SGB, species, Taxonomy)  ) -> NorthAmericaAss

EuropeAss %>% mutate(ID = paste0(SGB, "_", Phenotype ) ) %>% filter(ID %in%  mutate(AsiaAss, ID = paste0(SGB, "_", Phenotype ) )$ID )
EuropeAss %>% mutate(ID = paste0(SGB, "_", Phenotype ) ) %>% filter(ID %in%  mutate(NorthAmericaAss, ID = paste0(SGB, "_", Phenotype ) )$ID )



EuropeAss %>% group_by(Phenotype) %>% summarise(n())
AsiaAss %>% group_by(Phenotype) %>% summarise(n())
NorthAmericaAss %>% group_by(Phenotype) %>% summarise(n())


EuropeAss %>% group_by(SGB, species) %>% summarise(N = n()) %>% arrange(desc(N))
AsiaAss %>% group_by(SGB, species) %>% summarise(N = n()) %>% arrange(desc(N))
NorthAmericaAss %>% group_by(SGB, species) %>% summarise(N = n()) %>% arrange(desc(N))





###PCA####
library(tidyverse)
library(vegan)
library(phytools)
library(ape)

SGB = "t__SGB2318" #"t__SGB13977"
Tree = "/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB2318.TreeShrink.tre"
Metadata_location ="/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv" #"Phenotypes_merged.tsv"

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
pruned_tree <- keep.tip(Tree, Metadata$ID_anal )

#3.Removal of repeated samples
#Usually they cluster together. We could in principle leverage some random effects to accomodate multiple 
Unrepeated = Remove_repeated(Metadata, pruned_tree)
Metadata = Unrepeated[[1]] %>% ungroup() ; pruned_tree = Unrepeated[[2]]

#4. Make PCA
cophenetic.phylo( pruned_tree ) -> Phylo_dist
pcoa = cmdscale (Phylo_dist, eig = TRUE, k = dim(Metadata)[1]-1  )

Eigen_values = pcoa$eig / (pcoa$eig %>% sum())
Cumulative_variability_explained = cumsum(Eigen_values) 
PCs_check = which.max(Cumulative_variability_explained >= 0.9)

pcoa$points %>% as.data.frame() %>% rownames_to_column("ID_anal") %>% as_tibble() -> Components

Components %>% left_join(Metadata) %>% filter(! is.na(Age) )  -> Components

#Find countries with more than 1 study
Components %>% group_by(study_name, Country) %>% summarise(n()) %>% group_by(Country) %>% summarise(N = n())  %>% filter(N > 1) -> countries_multiple
#Find continents with more than 1 study
Components %>% group_by(study_name, Continent) %>% summarise(n()) %>% group_by(Continent) %>% summarise(N = n())  %>% filter(N > 1) -> continents_multiple
#Remove samples from studies with few samples
Components %>% group_by(study_name) %>% summarise(N = n()) %>%  filter(N < 10) -> studies_remove


ggplot(Components, aes(x=V1, y=V2, col=Continent)) + geom_point() + theme_bw()
tibble(N=seq(length(Cumulative_variability_explained)), Variability=Eigen_values) %>%  ggplot(aes(x=N,y=Variability)) + geom_bar(stat = "identity") + theme_bw()


All_results = tibble()

library(lmerTest)
Inverse_rank_normal = function(Measurement){
  qnorm((rank(Measurement,na.last="keep")-0.5)/sum(!is.na(Measurement)))
}

for (N in seq(PCs_check)){
  Component = paste0("V", N )
  Formula2 = paste0(Component, " ~ Age + Sex + BMI  + (1|Continent) + (1|Continent:Country) + (1|Continent:Country:study_name)")
  Components %>% filter(! study_name %in% studies_remove$study_name ) %>% filter(Country %in% countries_multiple$Country ) %>% filter(! (is.na(Age) | is.na(BMI) | is.na(Continent) | is.na(Sex) ) ) %>% lmer(Formula2, . ) -> Result
  confint(Result) -> CI
  Result %>% summary() -> Result
  Result$varcor %>% cbind(CI[1:3,] ) %>% as.data.frame() %>% rownames_to_column("Random_effect")  -> RE
  summary(Result)$varcor %>% as_tibble() %>% cbind(CI[1:dim(.)[1],]) %>% as_tibble() %>% select(-c(var1,var2)) %>% rename(Random_effect = grp)
  RE %>% mutate(Component= N, .before=1) %>% rbind(All_results, . ) -> All_results
  #Result$coefficients %>% as.data.frame() %>% rownames_to_column("Features") %>% as_tibble() %>% mutate(Component= Component) %>% rbind(All_results, . ) -> All_results
}




Formula2 = paste0(Component, " ~ Age + Sex + BMI + disease + (1|Continent/Country/study_name) ")
Components %>% filter(! study_name %in% studies_remove$study_name ) %>% filter(Country %in% countries_multiple$Country ) %>% lmer(Formula2, . ) %>% summary() -> Result



Components %>% filter(! study_name %in% studies_remove$study_name ) %>% filter(Country %in% countries_multiple$Country ) %>% MCMCglmm(V1 ~ Age +Sex + BMI + disease + (1|Contient/Country/study_name) , data = . )

prior1 = list(
  G=list(
    G1=list(V=1,nu = 0.001),
    G2=list(V=1,nu=0.001),
    G3=list(V=1,nu=0.001)),
  R=list(V=1,nu=0.001))

library(MCMCglmm)

All_results = tibble()

Get_info_from_MCMC = function(VCV){
  as_tibble(VCV) %>% summarise_all(median) -> info_median
  as_tibble(VCV) %>% summarise_all(mean) -> info_mean
  colnames(info_median) = paste0("median_", colnames(info_median) )
  colnames(info_mean) = paste0("mean_", colnames(info_mean) )
  posterior.mode(VCV) %>% t() %>% as.data.frame() %>% as_tibble() -> info1
  HPDinterval(VCV,prob = 0.95) %>% as.data.frame() %>% rownames_to_column('Variable')  %>%
    pivot_wider(names_from = Variable, values_from = c(lower, upper)) -> info2
  cbind(info1, info2) %>% cbind(info_median) %>% cbind(info_mean) %>% as_tibble() %>% mutate(Component = Component, .before=1 ) %>% return()
}

for (N in seq(PCs_check)){

  Components %>% filter(! study_name %in% studies_remove$study_name ) %>% filter(Country %in% countries_multiple$Country ) %>%
    select( - family ) %>% filter(!is.na(Age) ) %>% filter(!is.na(Sex) ) %>% filter(!is.na(BMI) )  -> Input
  
  Input$Continent_country <- paste0(Input$Continent,Input$Country)
  Input$Continent_country_study <- paste0(Input$Continent,Input$Country,Input$study_name)
    
   
   Input %>%  MCMCglmm::MCMCglmm(V2 ~ Age +Sex + BMI  ,
           random = ~Continent + Continent:Country+ Continent:Country:study_name, 
           prior = prior1,
           thin=50,
           data=. ,
           family = "gaussian",nitt = 60000,burnin=10000 ) -> result
   
   Input %>%  MCMCglmm::MCMCglmm(V2 ~ Age +Sex + BMI ,
                                 random = ~Continent + Continent_country+ Continent_country_study, 
                                 prior = prior1,
                                 thin=50,
                                 data=. ,
                                 family = "gaussian",nitt = 60000,burnin=10000 ) -> result2
   
   
  VCV2 = result$VCV
  VCV2 / rowSums(VCV2) -> VCV3
  #plot(VCV3)
  
  
  
  Get_info_from_MCMC(VCV2) %>%  rbind(All_results, . ) -> All_results
}

write_tsv(All_results, paste0("/mnt/project/Make_Associations/Association/Results/PCA_partion/", SGB, ".tsv") )






library(ape)
library(tidyverse)
library(ggtree)
library(phytools)
library(ggtreeExtra) #Allows to add metadata to circular plots
library(ggnewscale) #Allos for multiple fill scales within one plot
library(viridis) #color palette
source("/mnt/project/Make_Associations/Association/Functions.R")


Make_plot = function(Tree_file="/mnt/project/Symlink_phylogenies_with_Food_references/RAxML_bestTree.t__SGB5045.TreeShrink.tre", color_dataset=F){
#Tree = "/mnt/project/Symlink_phylogenies_with_Food_references/RAxML_bestTree.t__SGB7020.TreeShrink.tre"
SGB = str_split(Tree_file, "\\.")[[1]][2]
Plot_name = paste0(SGB,"_tree")

References = read_tsv("/mnt/project/Symlink_phylogenies_with_Food_references/References.txt", col_names = F)
print(References)



Metadata_location ="/mnt/project//Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv" #"Phenotypes_merged.tsv"

#1. Read. Tree
read.tree(Tree_file) -> Tree
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
References %>% filter(X1 %in% Tree$tip.label) -> References
#3. Read metadata
Metadata = read_tsv(Metadata_location) %>% filter(! study_name %in% c("ThomasAM_2019_c","500FG_FSK")  )
Metadata %>% group_by(ID_anal) %>% sample_n(1) -> Metadata
Metadata %>% filter(ID_anal %in% Tree$tip.label) -> Metadata
pruned_tree <- keep.tip(Tree, c(Metadata$ID_anal, References$X1  ))

#3.Removal of repeated samples
Remove_repeated = function(Metadata, pruned_tree){
  'Check which samples are repeated (either duplicates or longitudinally) and only picks one from the repeated'
  set.seed(89777)
  Metadata %>% group_by(subject_id) %>% summarise(N = n()) %>% arrange(desc(N)) %>% filter(N > 1 ) -> Replicates
  remove = tibble()
  for (i in Replicates$subject_id){
    Metadata %>% filter(subject_id == i) -> Filtered ; Samples = Filtered$ID_anal
    Keep = sample(Samples, 1)
    remove = c(remove, Samples[!Samples==Keep] ) 
  }
  remove = unlist(remove)
  Metadata %>% filter(! ID_anal %in% remove) -> Metadata
  pruned_tree %>% drop.tip(remove) -> pruned_tree
  return(list(Metadata, pruned_tree))
}
#Usually they cluster together. We could in principle leverage some random effects to accomodate multiple 
Unrepeated = Remove_repeated(Metadata, pruned_tree)
Metadata = Unrepeated[[1]] %>% ungroup() ; pruned_tree = Unrepeated[[2]]


#4. Make plots
Info = Metadata %>% select(ID_anal, study_name, Continent, Country, Age) 
Info %>% drop_na() -> Info
keep.tip(pruned_tree, c(Info$ID_anal, References$X1 ) ) -> pruned_tree
pruned_tree$tip.labe[!pruned_tree$tip.label %in% Info$ID_anal ] %>% as_tibble() %>% rename(ID_anal = value) %>%
  mutate(study_name = NA, Continent= NA, Country=NA, Age=NA, Is_reference=T ) -> Refs
Info %>% mutate(Is_reference=F) %>% rbind(Refs) -> Info

print("Making tree plot")
ggtree(pruned_tree, layout="fan", open.angle=15, size=0.1) %<+% Info -> p
p + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1 ) -> p2
p2 + scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) -> p2
p2 + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age), width=0.03,offset=0.1 ) -> p2
p2 + scale_fill_viridis_c(option = "viridis") -> p2
p2 + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Is_reference), width=0.03,offset=0.1 ) +  scale_fill_manual( values = c("grey", "#E83845") ) -> Plot
if (color_dataset == T){
  Plot + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Country), width=0.03,offset=0.1 )+  scale_fill_manual(values=c25)  + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=study_name), width=0.03,offset=0.1 )  +  scale_fill_manual(values=c25)
}
Output = str_replace(Tree_file, paste0("RAxML_bestTree.", SGB,".TreeShrink.tre"  ), paste0(SGB, ".pdf") )
print(Output)
ggsave(Output, Plot )
return(Plot)
}
#Make_plot(Tree="/mnt/project/Symlink_phylogenies_with_Food_references/RAxML_bestTree.t__SGB7142.TreeShrink.tre")

directory_path <- "/mnt/project/Symlink_phylogenies_with_Food_references/"
# List all files in the directory ending with ".tree"
tree_files <- list.files(path = directory_path, pattern = "\\.tre$", full.names = TRUE)
sapply(tree_files ,Make_plot)

###########################
##Proportion explained#####
###########################
library(tidyverse)
library(patchwork)


file_list <- list.files("/mnt/project/Make_Associations/Association/Results/PCA_partion/", pattern = "\\.tsv$", full.names = TRUE)
merged_data <- file_list %>% map_dfr(~read_tsv(.x) %>% mutate(File_Name = basename(.x)))
merged_data %>% mutate(SGB = str_replace(File_Name, "_noDisease_onlyHealthy.tsv", "") ) -> merged_data

#ICC calcualtion
merged_data %>% group_by(File_Name, Component) %>% mutate(IC = vcov / sum(vcov)) %>% ungroup() -> merged_data
#Multiply ICC by variance explained by component
merged_data %>% mutate(total_variance = IC * Variance_explained_PC ) -> merged_data

#Check proportions explained for "significant" study effect

merged_data %>%  filter(`2.5 %` > 0) %>% mutate(Random_effect = ifelse(Random_effect == "Continent:Country", "Country", ifelse(grepl("study", Random_effect), "Study", Random_effect ) ) ) %>% mutate(Random_effect = factor(Random_effect, levels =c("Residual", "Continent", "Country", "Study" ) ) ) -> For_plot
For_plot  %>%  ggplot(aes(y=IC, x=Random_effect )) +
  geom_half_boxplot(alpha=0.5,position = "identity", outlier.shape = NA) +  geom_half_violin(side = "r") + theme_bw() + coord_flip()  + ggtitle("Random effects variability explained (IC)") +  theme(plot.title = element_text(size = 7))   -> IC_distribution

For_plot %>% group_by(SGB, Random_effect) %>% summarise(total_variance_SGB = sum(total_variance)) %>%
  ggplot(aes(y=total_variance_SGB, x=Random_effect )) +
  geom_half_boxplot(alpha=0.5,position = "identity", outlier.shape = NA) +  geom_half_violin(side = "r") + theme_bw() + coord_flip() + ggtitle("Sum of IC*PC_proportion") +  theme(plot.title = element_text(size = 7)) -> Total_distribution

IC_distribution + Total_distribution -> final_plot
ggsave("/mnt/project/Make_Associations/Association/Results/PCA_partion/Variability_explained.pdf", final_plot)
write_tsv(merged_data, "/mnt/project/Make_Associations/Association/Results/PCA_partion/Merged_stats.tsv")



############Overlap Antropo associations
Asia_results = read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Asia_filtered.csv") %>% mutate(Key = paste0(Phenotype, "_", SGB ))
Europe_results = read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Europe_filtered.csv")%>% mutate(Key = paste0(Phenotype, "_", SGB ))
NA_results = read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica_filtered.csv")%>% mutate(Key = paste0(Phenotype, "_", SGB ))



Make_plot_comparison = function(Pheno){
  Europe_results %>% filter(Key %in% Asia_results$Key | Key %in% NA_results$Key )  %>% filter(Phenotype == Pheno)-> Check
  Asia_results %>% filter(Key %in% Check$Key ) %>% select(SGB, elpd_diff, phylo_mean ) %>% mutate(Continent = "Asia") %>% rbind(
    Europe_results %>% filter(Key %in% Check$Key ) %>% select(SGB, elpd_diff, phylo_mean ) %>% mutate(Continent = "Europe") ) %>% rbind(
      NA_results %>% filter(Key %in% Check$Key ) %>% select(SGB, elpd_diff, phylo_mean ) %>% mutate(Continent = "North_America") ) -> Check_stats
  Check_stats %>% ggplot(aes(x=SGB, y = elpd_diff, col = Continent )) + theme_bw() + geom_point() + coord_flip() + geom_hline(yintercept = -4, linetype=2) -> Plot1
  Check_stats %>% ggplot(aes(x=SGB, y = phylo_mean, col = Continent )) + theme_bw() + geom_point() + coord_flip() -> Plot2
  Plot1 + Plot2 -> Plot
  return(list(Plot, Check_stats))
}
Make_plot_comparison("Age") -> Res_Age
Make_plot_comparison("Sex")
Make_plot_comparison("BMI")

SGB_taxonomy = read_tsv("/mnt/project/Make_Associations/Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2", col_names = F )
SGB_taxonomy %>% rename(SGB=X1, Taxonomy=X2 ) %>% mutate(SGB = paste0("t__", SGB) ) -> SGB_taxonomy
Age_info = left_join(Res_Age[[2]],SGB_taxonomy)
Age_info <- Age_info %>%
  separate(Taxonomy, into = c("domain","phylum", "class", "order", "family", "genera", "species"), sep = "\\|", remove = FALSE)
Age_info %>% filter(Continent == "Europe")



Age_info %>% filter(SGB == "t__SGB17248" )



##Plot scatter tree size
read_tsv("/mnt/project/SGB_abundances/Result/Prevalence.tsv") -> Prevalences
Prevalences %>% filter(Continent == "All") -> Prevalences
str_split(Prevalences$SGB, "\\|") %>% lapply( function(x){ x[length(x)]} ) %>% unlist() ->SGB_names
Prevalences %>% rename( Taxonomy = SGB) %>% mutate(SGB = SGB_names) -> Prevalences
read_tsv("/mnt/project/Make_Associations/Association/Tree_sizes2.txt") -> tree_size
str_split(tree_size$Tree, "\\.") %>% lapply( function(x){ x[2]} ) %>% unlist() ->SGB_names
tree_size %>% mutate(SGB = SGB_names) %>% left_join(Prevalences) -> for_plot

str_split(for_plot$Taxonomy, "\\|") %>% lapply( function(x){ str_replace(x[2], "p__", "") } ) %>% unlist() -> phyla


lm( N_tips ~ Prevalence , for_plot ) -> linear_model 
exp_model <- nls(N_tips ~ a * exp(b * Prevalence), data = for_plot %>% drop_na(), start = c(a=74, b=6.8)  )  

for_plot %>% mutate(Phylum = phyla) %>% mutate(Phylum = ifelse(Phylum %in% c("Bacteroidetes", "Firmicutes","Proteobacteria", "Verrucomicrobia", "Lentisphaerae", "Actinobacteria"), Phylum, "Other" ) ) %>%
  ggplot(aes(x=Prevalence, y=N_tips )) + geom_point( aes(col=Phylum, alpha=ifelse(N_tips<300,0.2,0.6) ), size=1) + 
  geom_hline(yintercept = 300, linetype=2,  col= "black", linewidth=1 ) +
  theme_bw() + ylab("Number of samples\n in SGB phylogeny") + xlab("SGB global prevalence") +
  #geom_smooth(method="lm", formula= (y ~ log10(x)), se=FALSE, linetype = 1) +
  scale_color_manual(values = color_mapping) + scale_y_log10() + guides(alpha = FALSE) + theme(
    axis.text = element_text(size = 19), 
    axis.title = element_text(size = 19), 
    legend.text = element_text(size = 19),
    legend.title =  element_text(size = 19)
  ) -> Prevalence_vs_treeSize

for_plot %>% 
  mutate(Included = ifelse(N_tips >= 300, "Included", "Not_included")) %>% 
  group_by(Included) %>% 
  summarise(N. = n()) %>%
  ggplot(aes(x = "", y = N., fill = Included)) + 
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start = 0) + 
  theme_void() +  # Removes background and axes
  scale_fill_manual(values = c("Included" = "#E83845", "Not_included" = "grey")) + 
  labs(fill = "Tree included for analysis") +
  geom_text(aes(label = N.), 
            position = position_stack(vjust = 0.5),  # Positions labels in the middle of each slice
            color = "white", size = 5)  -> Piechart_inclusion


ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/Prevalence_vs_NPhylogeny.pdf", Prevalence_vs_treeSize, width = 6.5, height=4)
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/Prevalence_vs_NPhylogeny.png", Prevalence_vs_treeSize)  
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/Prevalence_vs_NPhylogeny.tiff", Prevalence_vs_treeSize)  
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/Inclusion_pie.pdf", Piechart_inclusion)  


  

color_mapping <- c(
  "Actinobacteria" = c25[1],
  "Firmicutes" = c25[2],
  "Bacteroidetes" = c25[3],
  "Proteobacteria" = c25[4],
  "Euryarchaeota" = c25[5],
  "Synergistetes" = c25[7],
  "Lentisphaerae" = c25[8],
  "Kiritimatiellaeota" = c25[9],
  "Tenericutes" = c25[10],
  "Elusimicrobia" = c25[11],
  "Crenarchaeota" = c25[13],
  "Spirochaetes" = c25[14],
  "Fusobacteria" = c25[15],
  "Verrucomicrobia" = c25[15],
  "Planctomycetes" = c25[17],
  "Bacteria_unclassified" = c25[18],
  "Candidatus_Melainabacteria" = c25[19],
  "Candidatus_Thermoplasmatota" = c25[20],
  "Candidatus_Saccharibacteria" = c25[21],
  "Other" = c25[22]
)
c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")



#####
library(ggrepel)
read_tsv("/mnt/project/Strain_sharedness/Tranmission_rates.tsv") -> df_transmission
Phylo = read_tsv("/mnt/project/Make_Associations/Association/Results/Geography_cor/Association_results_merged.tsv")
df_transmission %>% mutate(Relationship = ifelse(motherbaby==T, "MotherBaby", ifelse(Family==T, "Family", ifelse(Intercountry==T,"Intercountry", 
  ifelse(  Intracountry_diffStudy == T, "Intracountry_diffStudy", ifelse(Intracountry_sameStudy == T, "Intracountry_sameStudy" , NA ) ) ) ) ) ) -> df_transmission
df_transmission %>% select(Relationship, Shared, SGB) %>% spread(Relationship, Shared) -> df_transmission

print("===Missing SGBS=====")
Phylo %>% filter(!SGB %in% df_transmission$SGB) %>% print()
left_join(df_transmission, Phylo) -> df_transmission

#Mark E. rectale (SGB4933_group), R. torques (SGB4563_group), other interesting examples? 
#In ordre: Troques (clear caldes),  rectale (previously seen in Nicola), aerofaciens (highlighted in Suzuki), Ventrisium (large effect in more than 3K), Lachnospiraceae (top oen with >6K), Alistipes_onderdonkii (low one with over 10K samples)
SGB_geography_highlight  = c("SGB4563_group", "SGB4933_group", "SGB14546_group", "SGB5045", "SGB4910", "SGB2303" ) 

df_transmission %>% lm(Rho ~ Intracountry_diffStudy , . ) %>% summary()
df_transmission %>% filter(FDR<0.05) %>% mutate(SGB = str_replace(SGB, "t__", "") ) %>% mutate(Bug = ifelse(SGB == "SGB4933_group", "E.rectale", ifelse(SGB == "SGB4563_group", "R.torques", ifelse(SGB=="SGB14546_group", "C.aerofaciens", 
                                                                                       ifelse(SGB=="SGB5045", "E.ventrisium", ifelse(SGB=="SGB4910", "unknown_Lachnospiraceae" , ifelse(SGB=="SGB2303", "A.onderdonkii", NA) ) ) )  ) ) ) %>%
  ggplot(aes(x=Intracountry_diffStudy, y= Rho ) ) +  geom_point(aes(col= SGB %in% SGB_geography_highlight)) + theme_bw() +
  geom_smooth(method = "lm") + geom_label_repel(aes(label = Bug), data = . %>% filter(SGB %in% SGB_geography_highlight)) +
  xlab("Strain sharing proportion between\nsamples from same country but different study") + ylab("Geography effect (rho)") + scale_color_manual(values= c("black", "red") ) + theme(legend.position = "none", text = element_text(size = 16)) -> Plot_sharingVSGeography
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/StrainSharinv_vs_GeogrpahyRho_sign.tiff" ,Plot_sharingVSGeography)
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/StrainSharinv_vs_GeogrpahyRho_sign.pdf" ,Plot_sharingVSGeography)

cor.test(df_transmission$Intracountry_diffStudy, df_transmission$Rho, method = "spearman", na.action="omit" )



###Plot trees bilirubin
Make_plot_bilirubin = function(Tree_file="/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB5045.TreeShrink.tre", pheno= "bilirubin", with_disease=T, Cont = "All", Model_name="No" ){
  #Tree = "/mnt/project/Symlink_phylogenies_with_Food_references/RAxML_bestTree.t__SGB7020.TreeShrink.tre"
  SGB = str_split(Tree_file, "\\.")[[1]][2]
  Plot_name = paste0(SGB,"_tree")
  
  Metadata_location ="/mnt/project//Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv" #"Phenotypes_merged.tsv"
  
  #1. Read. Tree
  read.tree(Tree_file) -> Tree
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
  References %>% filter(X1 %in% Tree$tip.label) -> References
  #3. Read metadata
  Metadata = read_tsv(Metadata_location) %>% filter(! study_name %in% c("ThomasAM_2019_c","500FG_FSK")  )
  Metadata %>% group_by(ID_anal) %>% sample_n(1) -> Metadata
  Metadata %>% filter(ID_anal %in% Tree$tip.label) -> Metadata
  pruned_tree <- keep.tip(Tree, c(Metadata$ID_anal, References$X1  ))
  
  #3.Removal of repeated samples
  Remove_repeated = function(Metadata, pruned_tree){
    'Check which samples are repeated (either duplicates or longitudinally) and only picks one from the repeated'
    set.seed(89777)
    Metadata %>% group_by(subject_id) %>% summarise(N = n()) %>% arrange(desc(N)) %>% filter(N > 1 ) -> Replicates
    remove = tibble()
    for (i in Replicates$subject_id){
      Metadata %>% filter(subject_id == i) -> Filtered ; Samples = Filtered$ID_anal
      Keep = sample(Samples, 1)
      remove = c(remove, Samples[!Samples==Keep] ) 
    }
    remove = unlist(remove)
    Metadata %>% filter(! ID_anal %in% remove) -> Metadata
    pruned_tree %>% drop.tip(remove) -> pruned_tree
    return(list(Metadata, pruned_tree))
  }
  #Usually they cluster together. We could in principle leverage some random effects to accomodate multiple 
  Unrepeated = Remove_repeated(Metadata, pruned_tree)
  Metadata = Unrepeated[[1]] %>% ungroup() ; pruned_tree = Unrepeated[[2]]
  if (Cont != "All"){ Metadata %>% filter(Continent == Cont) -> Metadata }
  
  
  #4. Make plots
  if (with_disease == T){
   Info = Metadata %>% select( c("ID_anal", "study_name", "Continent", "Country", "Age", pheno, "disease") ) 
  } else {  Info = Metadata %>% select( c("ID_anal", "study_name", "Continent", "Country", "Age", pheno, "BMI") )   }
  Info %>% drop_na() -> Info
  keep.tip(pruned_tree, c(Info$ID_anal ) ) -> pruned_tree
  
  #Add phylogenetic signal if needed
  if (Model_name != "No"){
    Model = read_rds(Model_name)
    Model$model_input -> Info_model
    Model$pglmm_fit$summary() -> Fit
    Fit %>% filter( grepl("phylo", variable))  %>% filter(! grepl("std", variable)) %>% mutate(sample_id = factor(c("-", as.character(Info$sample_id)))) %>%
    select(sample_id, median) %>% rename(phylo_effect_median = median) -> Phylo
    left_join(Info_model, Phylo) %>% mutate(ID_anal = as.character(sample_id)) -> Info_model
    left_join(Info ,Info_model) -> Info
  }
  
  
  print("Making tree plot")
  ggtree(pruned_tree, layout="fan", open.angle=15, size=0.1) %<+% Info -> p
  if (Cont == "All"){
    p + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1 ) -> p2
    p2 + scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) -> p2
  } else { p + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Country), width=0.03,offset=0.1 ) -> p2  }
  p2 + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age), width=0.03,offset=0.1 ) -> p2
  p2 + scale_fill_viridis_c(option = "viridis") -> p2
  if (with_disease == T){
  color_mapping_di <- c(
      "T2D" = "#A6CEE3",
      "T2D;hypertension" = "#1F78B4",
      "ascites;cirrhosis" = "#B2DF8A",
      "ascites;cirrhosis;hepatitis" ="#33A02C",
      "ascites;cirrhosis;hepatitis;schistosoma" = "#33A02C",
      "ascites;cirrhosis;hepatitis;wilson" = "#33A02C",
      "cirrhosis" = "#B2DF8A",
      "cirrhosis;hepatitis " = "#33A02C",
      "healthy" ="white",
      "hepatitis" = "#FDBF6F",
      "hypertension" ="#FB9A99"
    )
    p2 + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=disease), width=0.03,offset=0.1 ) +  scale_fill_manual( values = color_mapping_di ) -> Plot
  } else { p2 + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=BMI), width=0.03,offset=0.1 ) +  scale_fill_viridis_c(option = "mako") -> Plot  }
  Plot + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Country), width=0.03,offset=0.1 )+  scale_fill_manual(values=c25)  + new_scale_fill() -> Plot
  Plot + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=study_name), width=0.03,offset=0.1 )  +  scale_fill_manual(values=c25) -> Plot
  if (Model_name != "No"){
    Plot + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") -> Plot
  }  
  Plot +  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes_string(fill=pheno), width=0.03,offset=0.1 ) 

}
Make_plot_bilirubin("/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB1903_group.TreeShrink.tre")

#blood pressure
Make_plot_bilirubin("/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB3926_group.TreeShrink.tre", 
pheno = "SBP", with_disease = F, Cont = "Europe", 
Model_name = "/mnt/project/Make_Associations/Association/Results/Continent_stratified/Age,BMI_Europe/Models/t__SGB3926_group/SBP/Model.rds"  )



####Get per phnotypic category how many phenotypes

read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypic_categories.tsv") -> Pheno
To_keep = read_tsv("/mnt/project/Make_Associations/Phenotypes/Remove_variables.tsv")

Pheno %>% filter(! Phenos %in% To_keep$Remove) %>% group_by(Category) %>% summarise(n())



###############################################################
##############Correlation network only in DAG3#################
###############################################################

#We can potentially replicate in other studies...

#1. Get DAG3 samples
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") %>% filter(study_name == "DAG3") %>% select(ID_anal, study_name) -> Samples_dag3
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") %>% filter(study_name == "MetaCardis_2020_a") %>% select(ID_anal, study_name) -> Samples_metacardis
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") %>% filter(study_name == "AsnicarF_2021") %>% select(ID_anal, study_name) -> Samples_asnicar
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") %>% filter(study_name == "ZeeviD_2015") %>% select(ID_anal, study_name) -> Samples_zheevi



#2. Functions
Keep_common_nodes = function(Tree1, Tree2){
  
}
Prepare_tree = function(Tree, ID_list){
  read.tree(Tree) -> Tree
  if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }
  Tree = midpoint.root(Tree)
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
  
  ID_list %>% filter(ID_anal %in% Tree$tip.label) -> IDs
  Tree = keep.tip(Tree , IDs$ID_anal)
  
  return(Tree)
}

Keep_common_nodes = function(Tree1, Tree2){
  Names = c(Tree1$tip.label, Tree2$tip.label)
  Names = Names[duplicated( Names )] %>% unique()
  Tree1 = keep.tip(Tree1 , Names)
  Tree2 = keep.tip(Tree2 , Names)
  return(list(Tree1, Tree2))
}

Get_SGB = function(Path){
  str_split( basename(Path) , "\\." )[[1]][2]
}

Correlate_trees = function(Tree1n, Tree2n, ID_list=Samples_dag3){
  Prepare_tree(Tree1n,ID_list )  -> Tree1
  Prepare_tree(Tree2n,ID_list )  -> Tree2
  if (is.null(Tree2) | is.null(Tree1)){ return(tibble(Tree1 = Tree1n , Tree2 = Tree2n , Common_samples=0, Rho= NA, P=NA )) }

  Trimmed = Keep_common_nodes(Tree1, Tree2)
  Tree1 = Trimmed[[1]]
  Tree2 = Trimmed[[2]]
  N = length(Tree1$tip.label)
  
  if (N < 20 ) { return(tibble(Tree1 = Tree1n , Tree2 = Tree2n , Common_samples=N, Rho= NA, P=NA )) } 
  
  cophenetic.phylo(Tree1) -> D1
  cophenetic.phylo(Tree2) -> D2
  D1[rownames(D2), colnames(D2)] -> D1
  
  vegan::mantel( xdis= D1, ydis = D2 ) -> res
  Rho = res$statistic
  P = res$signif
  return( tibble(Tree1 = Get_SGB(Tree1n) , Tree2 = Get_SGB(Tree2n) , Common_samples=N, Rho= Rho, P=P ) )
}

Tree1n = "/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB9281.TreeShrink.tre"
Tree2n = "/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB9273.TreeShrink.tre"
Correlate_trees(Tree1n, Tree2n, Samples_dag3)

read_tsv("/mnt/project/Make_Associations/Association/Tree_sizes2.txt") %>% filter(N_tips > 1000) %>% arrange(N_tips) -> For_analysis
Results_corr = tibble()
for (N1 in seq(length(For_analysis$Tree)) ){
  if (N1 == length(For_analysis$Tree)){ next }
  Tree1n = paste0("/mnt/project/Symlink_phylo/", For_analysis$Tree[N1] )
  for (N2 in seq(N1+1, length(For_analysis$Tree) ) ){
    Tree2n = paste0("/mnt/project/Symlink_phylo/", For_analysis$Tree[N2] )  
    print(paste0(Tree1n, " ", Tree2n) )
    Results_corr = rbind(Results_corr, Correlate_trees(Tree1n, Tree2n, Samples_asnicar) )
  }
}
write_tsv(Results_corr, "/mnt/project/Make_Associations/Association/Results/Correlation_trees/Summaries_asnicar.tsv")



Abundance %>% select(-ID) %>% apply(2, function(x){ sum(x!=0)/length(x) }  ) %>% tibble(SGB = names(.), Prevalence = . ) -> Prevalence
Prevalence %>% filter(!Prevalence < 0.1 ) -> Keep_SGB




##############################
#Analyze correlation trees####
##############################

c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
color_mapping <- c(
  "Actinobacteria" = c25[1],
  "Firmicutes" = c25[2],
  "Bacteroidetes" = c25[3],
  "Proteobacteria" = c25[4],
  "Euryarchaeota" = c25[5],
  "Synergistetes" = c25[7],
  "Lentisphaerae" = c25[8],
  "Kiritimatiellaeota" = c25[9],
  "Tenericutes" = c25[10],
  "Elusimicrobia" = c25[11],
  "Crenarchaeota" = c25[13],
  "Spirochaetes" = c25[14],
  "Fusobacteria" = c25[15],
  "Verrucomicrobia" = c25[15],
  "Planctomycetes" = c25[17],
  "Bacteria_unclassified" = c25[18],
  "Candidatus_Melainabacteria" = c25[19],
  "Candidatus_Thermoplasmatota" = c25[20],
  "Candidatus_Saccharibacteria" = c25[21],
  "Other" = c25[22]
)
tibble(C = color_mapping, Taxa = names(color_mapping) ) -> color_mapping


SGB_taxonomy = read_tsv("/mnt/project/Make_Associations/Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2", col_names = F ) %>%   separate(X2, into = c("domain","phylum", "class", "order", "family", "genera", "species"), sep = "\\|", remove = FALSE) %>% mutate(SGB = paste0("t__", X1))
directory_path = "/mnt/project/Make_Associations/Association/Results/Correlation_trees/"


Build_network = function(dataset="dag3"){

tsv_files <- list.files(directory_path, pattern = paste0("Summaries_",dataset,"_.*\\.tsv"), full.names = TRUE)
Correlation_df <- lapply(tsv_files, read_tsv) %>% bind_rows()
Correlation_df %>%  mutate(ID = ifelse(Tree1 < Tree2, paste(Tree1, Tree2, sep = "_"), paste(Tree2, Tree1, sep = "_"))) -> Correlation_df
Correlation_df %>% drop_na() %>% distinct(ID, .keep_all = TRUE)  -> Correlation_df2


Correlation_df2 %>% filter(Common_samples >= 200) %>% arrange(desc(abs(Rho))) %>% mutate(FDR = p.adjust(P, "fdr"))  -> Correlation_df_analysis
  
library(igraph)  
Correlation_df_analysis %>% filter(P==0.001, abs(Rho) > 0.1  ) -> Correlated_trees
graph_from_data_frame( Correlated_trees, directed = F ) -> Graph_cor

# Set edge widths based on the "Rho" attribute
edge_weights <- E(Graph_cor)$Rho

names(V(Graph_cor)) -> N
SGB_taxonomy[ match(N, SGB_taxonomy$SGB), ]  -> Taxonomies_color
Taxonomies_color %>% select(phylum)
color_mapping %>% mutate(phylum = paste0("p__", Taxa) ) %>% left_join(Taxonomies_color,.) %>% mutate(C=ifelse(is.na(C), "yellow4"  ,C) ) -> Taxonomies_color2

Taxonomies_color2 %>% group_by(family) %>% summarise(n()) %>% dim() -> N_till
Taxonomies_color2 %>% group_by(family) %>% summarise(n()) %>% select(family) %>% mutate(C2 = c25[1:N_till[1]] ) %>% left_join(Taxonomies_color2, . ) -> Taxonomies_color2
Taxonomies_color2$species = Taxonomies_color2$species %>% sapply(function(x){str_split(x, ",")[[1]][1] } )


#Correlation_df_analysis

# Plot the graph with variable edge widths
set.seed(2342)
plot(
  Graph_cor,
  edge.width = edge_weights * 10,  # Adjust the multiplier to control edge thickness
  vertex.size = 5,              # Optional: set the size of the vertices
  layout = layout_nicely(Graph_cor),  # Optional: arrange the layout
  vertex.color = Taxonomies_color2$C2,
  vertex.label.color = Taxonomies_color2$C,
  vertex.label = Taxonomies_color2$species,
  vertex.label.cex = 0.7
)

degree_centrality <- degree(Graph_cor)
Centrality_df = tibble( SGB = names(degree_centrality), DC = degree_centrality )
left_join(Centrality_df, Taxonomies_color2) %>% arrange(desc(DC)) %>% select(SGB, DC, species, genera, family) %>% print()

return()
#1 f__Acidaminococcaceae                   2 "dodgerblue2"
#2 f__Alphaproteobacteria_unclassified     1 "#E31A1C"    
#3 f__Bacilli_unclassified                 1 "green4"     
#4 f__Bacteroidales_unclassified           2 "#6A3D9A"    
#5 f__Clostridia_unclassified              9 "#FF7F00"    
#6 f__Clostridiaceae                       1 "black"     
#7 f__Clostridiales_unclassified           1 "gold1"      
#8 f__Lachnospiraceae                      5 "skyblue2"  
#9 f__Oscillospiraceae                     1 "#FB9A99"   
#10 f__Pirellulaceae                        1 "palegreen2" 
#11 f__Prevotellaceae                       5 "#CAB2D6"    
#12 f__Rikenellaceae                        2 "#FDBF6F"    
#13 f__Ruminococcaceae                     20 "gray70"     

}

Build_network("MetaCardis_2020_a")
Build_network("ZeeviD_2015")
Build_network("AsnicarF_2021")


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
#write_tsv(Results_enrichment_p, "/mnt/project/Make_Associations/Association/Results/Summaries/Enrichment_taxonomiclevel_associations.tsv")
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



################################
#Summarise with plots###########
################################
library(UpSetR)

#1. By phenotype
All_associations %>% group_by(Phenotype, Continent) %>% summarise(N = n()) %>% spread(Continent, N) %>% ungroup() -> N_phenotypes
N_phenotypes[is.na(N_phenotypes)] = 0
N_phenotypes %>% left_join(Categories) %>% select(-Phenotype) %>% group_by(Category) %>%
  summarise_all(sum) %>% gather(Continent, N_associations,Africa:South_America) %>% mutate(Category = ifelse(Category == "Lifestyle_and_exposome", "Lifestyle &\nexposome", ifelse(Category == "Medication_and_Supplements", "Medication &\n supplements",  Category) ) ) %>% drop_na() %>%
  ggplot(aes(x=reorder(Category, N_associations), y=N_associations, fill=Continent )) + geom_bar(stat="identity") + theme_bw() +
  scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) + coord_flip() +
  xlab("Phenotypic category") + theme(
    text = element_text(size = 19), # Adjust the text size as needed
    axis.text.x = element_text(size = 19), # Text size for x-axis labels
    axis.text.y = element_text(size = 19), # Text size for y-axis labels
    legend.text = element_text(size = 19) # Text size for legend
  ) + ylab("Association number") -> Plot1
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/BarPlot_AssociationN.tiff",Plot1)  
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/BarPlot_AssociationN.pdf",Plot1,width = 8, height=5.5)  

  
#2. By common associations
All_associations %>% select(ID, Continent) %>% mutate(Present = 1) %>% spread(Continent, Present ) -> All_associations_wide
All_associations_wide[is.na(All_associations_wide)] = 0
upset(All_associations_wide %>% as.data.frame(), sets.x.label = "SGB tree-Phenotype\nsupported assocations", 
      mainbar.y.label = "Number intersected supported associatons", text.scale = 2, nsets = 6, nintersects = 15 ) -> Plot2
pdf(file="/mnt/project/Make_Associations/Association/Results/Plots/UpSet_AssociationN.pdf") # or other device
Plot2
dev.off()


All_associations_wide %>% filter(Asia == 1 & Europe ==1) %>%
  separate(ID, into = c("Phenotype", "SGB"), sep = "-", extra = "merge") %>% group_by(Phenotype) %>% summarise(n())





############################
###Phenotype association####
############################

#Comparison abundance and tree signal
Prefix ="/mnt/project/"
Prefix ="/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/" #"/mnt/project/"


library(tidyverse)
library(MCMCglmm)
library(ape)
library(phytools)
source("Functions.R")


read_tsv( paste0(Prefix,"Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv")) %>% filter(study_name == "DAG3") -> Gacesa_phenos
For_analysis = c("Age", "BMI", "HDL", "LDL", "Appendectomy", "PPI", "NoDisease", "PM2.5")
Gacesa_phenos %>% select( c("ID_anal", For_analysis) ) -> Gacesa_phenos


#Abundance
read_rds(paste0(Prefix,"SGB_abundances/Data/Merged_Groningen_SGB.rds")) %>% filter(ID %in% Gacesa_phenos$ID_anal) %>% rename(ID_anal = ID) -> Abundance
Abundance %>% distinct(ID_anal, .keep_all = T) -> Abundance

Abundance %>% select(-ID_anal) %>% apply(2, function(x){ sum(x!=0)/length(x) }  ) %>% tibble(SGB = names(.), Prevalence = . ) -> Prevalence
Prevalence %>% filter(!Prevalence < 0.1 ) -> Keep_SGB
Abundance %>% select( c("ID_anal", Keep_SGB$SGB) ) -> Abundance
#Apply CLR
print("CLR-transformation")
Geom_mean <- function(x){
  exp(mean(log(x)))
}
CLR <- function(D){
  log2(D / Geom_mean(D))
}
Pseudo = min(select(Abundance, -ID_anal)[select(Abundance, -ID_anal) > 0 ] ) / 2
Abundance %>% select(-ID_anal) %>% apply(1, function(x){ CLR(x + Pseudo) } ) %>% t() %>% as_tibble() -> CLR_transformed
colnames(CLR_transformed) = colnames(select(Abundance, -ID_anal))
CLR_transformed %>% select(-UNCLASSIFIED) %>% mutate(ID_anal = Abundance$ID_anal ) -> CLR_transformed
CLR_transformed %>% select(-ID_anal) %>% colnames() %>% sapply( function(x){str_split(x, "\\|")[[1]] -> y ; y[length(y)]  }  ) -> Short_names
colnames(CLR_transformed) = c(Short_names, "ID_anal")


Prevalence$SGB %>% sapply( function(x){str_split(x, "\\|")[[1]] -> y ; y[length(y)]  }  ) -> Prevalence$SGB
Prevalence %>% arrange(desc(Prevalence)) %>% filter(SGB != "UNCLASSIFIED") %>% head(n=20) -> Most_prevalent

Process_tree = function(Tree){
  read.tree(Tree) -> Tree
  if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }
  Tree = midpoint.root(Tree)
  #2. Rename labels tree
  print("Processing tree")
  Tree$tip.label %>% sapply(. , Clean_names) -> New_names
  Tree$tip.label = as.vector(New_names)
  #if ( length(Tree$tip.label) < 100 ) { print("Tree is too small") ; q() }
  return(Tree)
  
}


prior1<-list(G=list(G1=list(V=1,nu=0.01)),R=list(V=1,nu=0.01))
prior0<-list(R=list(V=1,nu=0.01))


data4analysis = left_join(Gacesa_phenos, CLR_transformed)  %>% drop_na()                         
for ( Sp in as.vector(Most_prevalent$SGB) ) {
  #find tree
  Tree = paste0(paste0(Prefix,"Symlink_phylo/IQtree.",Sp, ".TreeShrink.tre"))
  if (! file.exists(Tree)){ Tree =paste0(paste0(Prefix,"Symlink_phylo/RAxML_bestTree.",Sp, ".TreeShrink.tre"))  }
  Process_tree(Tree) -> Tree
  data4analysis %>% filter(ID_anal %in% Tree$tip.label) -> data4analysis_f
  keep.tip(Tree, data4analysis_f$ID_anal) -> Tree
  #Ainv.1 = inverseA(Tree,nodes = "TIPS",scale = F)$Ainv
  data4analysis_f =  as.data.frame(data4analysis_f)
  for ( Pheno in  For_analysis){
    COV = c(Sp,"Age")
    COV = COV[! COV == Pheno ]
    #Formula = as.formula(paste0(Pheno, " ~ ", Sp, paste( c("",COV), collapse"+")  ))

    
    anpan_pglmm(data4analysis_f %>% rename(sample_id = ID_anal) , Tree, outcome = Pheno, covariates = COV, omit_na = T, family="gaussian", save_object =F,parallel_chains = 4, loo_comparison=T, show_plot_cor_mat=F, show_plot_tree=F, show_post=F ) -> result
    
    #Is the model better?
    result$loo$comparison %>% as_tibble() -> LOO
    LOO[1,1:2] -> LOO
    #How does the estimate change
    Fit1 = result$pglmm_fit$summary() %>% filter(variable == "beta[1]") %>% mutate(Model ="phylogeny") %>% mutate(Phenotype = Pheno, SGB=Sp)
    Fit0 = result$base_fit$summary() %>% filter(variable == "beta[1]") %>% mutate(Model = "base") %>% mutate(Phenotype = Pheno, SGB=Sp)
    rbind(Fit0, Fit1) %>% mutate(elpd_diff = LOO$elpd_diff, se_diff = LOO$se_diff) -> Model_results
    
    
    
#DIC will be in summary(MCMCglmm.object)
#posterior effect sizes of fixed effects:  MCMCglmm.object$Sol
  }
  
}



glm(Formula , data = data4analysis_f, family = "gaussian") -> Model0 #, cov_ranef = list(sample_id = Tree), s2.init = c(BM$ss)^2 ) -> CM
phyr::pglmm( Formula2, data = data4analysis_f, family = "gaussian", REML = T, cov_ranef = list(ID_anal = Tree ) ) -> CM 

#anpan_pglmm(data4analysis_f %>% rename(sample_id=ID_anal ) , Tree, outcome = Pheno, covariates = Sp, omit_na = T, family="gaussian", save_object =F,parallel_chains = 4, loo_comparison=F ) -> result



#####################
#Bifido animalis#####
#####################

#read.tree("/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB17278.TreeShrink.tre") -> Bifido_animalis
#read.tree("/mnt/project/Symlink_phylogenies_with_Food_references/RAxML_bestTree.t__SGB17278.TreeShrink.tre") -> Bifido_animalis2
read.tree("/mnt/project/Phylogenies/t__SGB17278_withFood/RAxML_bestTree.t__SGB17278.StrainPhlAn4.tre") -> Bifido_animalis3

read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> Phenotypes
Clean_names =  function(Name){
  'Samples from 500FG FSK need everything after _ to be removed. If _metaphlan4 is attached at the end of the name, remove'
  if ( grepl("HV", Name)) {
    Name = str_split(Name, "_")[[1]][1]
  }   
  Name = str_replace(Name, "_metaphlan4", "") 
  return(Name)
}
Bifido_animalis3$tip.label %>% sapply(. , Clean_names) %>% as.vector() -> Bifido_animalis3$tip.label
Phenotypes %>% filter(ID_anal %in% Bifido_animalis3$tip.label  )  -> Phenotypes_ba

Bifido_animalis3$tip.label[!Bifido_animalis3$tip.label %in% Phenotypes_ba$ID_anal] -> Missing

Phenotypes_ba %>% select(ID_anal, Country, Continent, study_name ) ->Phenotypes_ba
Phenotypes_ba %>% mutate(Reference_food = F, Mice=F ) %>% rbind( tibble(ID_anal =Missing   , Country=rep(NA, length(Missing)) , Continent=rep(NA, length(Missing)), study_name=rep(NA, length(Missing)), Reference_food=c(T, rep(NA, length(Missing)-1)), Mice=c(F,rep(T, length(Missing)-1) )  ) ) -> Phenotypes_ba
#M1231812230     ER_001__ACT__bin.2      ACT     ER      MAG     none    none    SGB17278        GGB10640        FGB3150 60.22   0.0     0.0     http://cmprod1.cibio.unitn.it/databases/MetaRefSGB/resources/Jan21/sequences/00000631/M1231812230.fna.bz2
#ER      none    Environmental   Food production Dairy products  Fermented dairy products        Yoghurt and dietary supplement  none    none    none    none    none    none


#Distances
cophenetic.phylo(Bifido_animalis3) -> Distances
dist_df <- as.data.frame(as.table(Distances))
colnames(dist_df) <- c("rowname", "colname", "distance")
as_tibble(dist_df) %>% filter(!rowname==colname) -> dist_df
#Add info
left_join(dist_df, Phenotypes_ba, by=c("rowname" = "ID_anal") ) %>% left_join(., Phenotypes_ba, by=c("colname" = "ID_anal"),  suffix=c("_1", "_2" ) ) -> dist_df
dist_df %>% mutate(Same_continent = ifelse(Continent_1 == Continent_2, T, F), Same_country = ifelse(Country_1 == Country_2, T, F), Same_study = ifelse(study_name_1 == study_name_2, T, F)  ) -> dist_df

dist_df %>%  ggplot( aes(x= as.vector(distance), fill=Same_study ) ) + geom_histogram() + theme_bw() + scale_x_log10() + xlab("Phylogenetic distance")

filter(dist_df %>% filter(distance< 1e-3) )  -> Cluster

Phenotypes_ba %>% mutate(Clade = ifelse(ID_anal %in% c(Cluster$rowname, Cluster$colname), T, F  ), ) -> Phenotypes_ba


Phenotypes_ba %>% filter(!study_name %in% c("DAG3", "LLD")) -> Phenotypes_ba2
keep.tip(Bifido_animalis3, Phenotypes_ba2$ID_anal) -> Bifido_animalis3
ggtree(Bifido_animalis3, open.angle=15, size=0.1)  %<+% (Phenotypes_ba2) -> p
p + geom_tippoint( mapping = aes(col=Reference_food)) +  scale_color_manual( values = c( "FALSE"="grey", "TRUE" = "#E83845") ) -> p
p + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1 ) -> p #+ scale_fill_manual( values = c( "FALSE"="grey", "TRUE" = "#E83845") )
p + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=study_name), width=0.03,offset=0.1 ) #+ scale_fill_manual( values = c( "FALSE"="grey", "TRUE" = "#E83845") )


dist_df %>% filter( (Mice_1 == T | Mice_2 == T) & ! (Mice_1 == Mice_2)   ) %>% mutate(Distance = "Human_to_Mouse" ) -> Distances_to_mice
dist_df %>% filter( (Mice_1 == F & Mice_2 == F) &  (Reference_food_1 == T |  Reference_food_2 == T) &  !(Reference_food_1 == Reference_food_2 )   ) %>% mutate(Distance = "Human_to_Food") -> Distances_to_food
dist_df %>% filter( (Mice_1 == F & Mice_2 == F) & !( Reference_food_1 == T |  Reference_food_2 == T )   ) %>% mutate(Distance = "Human_to_Human") -> Distances_human
dist_df %>% filter( (Mice_1 == T & Mice_2 == T)   ) %>%  mutate(Distance = "Human_to_Mouse") -> Distances_mice

rbind(rbind(rbind(Distances_to_mice, Distances_to_food), Distances_human), Distances_mice) %>% mutate( Relationship = ifelse( (Continent_1 != Continent_2 & Distance == "Human_to_Human" ) , "Different continent", 
                                                                                                                                ifelse( (Country_1 != Country_2 & Distance == "Human_to_Human"), "Same continent; different country",  
                                                                                                                                ifelse((study_name_1 != study_name_2 & Distance == "Human_to_Human" ), "Same country; different study",                                                                                 
                                                                                                                                ifelse((study_name_1 == study_name_2 & Distance == "Human_to_Human" ), "Same study" , NA))))) -> Distances_to_plot
Distances_to_plot %>% mutate( Distance =  ifelse(Mice_1 == T & Mice_2 == T, "Mouse-Mouse", Distance) ) %>% filter(! Distance == "Mouse-Mouse" ) %>% 
  ggplot(aes(x=Distance, y =distance, fill=Relationship )) + geom_boxplot() + theme_bw() + scale_y_continuous(trans='log10') + ylab("Phylogenetic distance") + xlab("Distance between") + scale_fill_manual( values= c("Different continent" = "#E31A1C", "Same continent; different country"= "#FDBF6F", "Same country; different study"= "steelblue4", "Same study"= "darkgreen"  )  ) + coord_flip() +
  theme(
    text = element_text(size = 16),  # Set the text size for all elements
    plot.title = element_text(size = 20, face = "bold"),  # Set title text size
    axis.title = element_text(size = 18, face = "bold"),  # Set axis title text size
    axis.text = element_text(size = 14),  # Set axis text size
    axis.ticks = element_line(size = 1.5),  # Set axis tick size
    legend.text = element_text(size = 14),  # Set legend text size
    legend.title = element_text(size = 16, face = "bold")  # Set legend title text size
  ) -> Plot_distances_banimalis
ggsave("/mnt/project/Make_Associations/Git/PhylogeneticAssociation/Analyses_revision/Figures/SF1_D.tiff", Plot_distances_banimalis , width = 8 , height= 6, units = "in" )

Distances_to_plot %>% filter(Distance == "Human_to_Human") %>% select(distance, Relationship) -> Test_human_diff
Test_human_diff %>% aov(distance ~ Relationship ,.) -> R
summary(R)[[1]][4][[1]][1] -> F_value
#Permutations
N_permutation = 1000
Perm_Fs = c()
set.seed(22311)
for (i in seq(1,N_permutation)){
  Test_human_diff_p = Test_human_diff 
  Perm = sample(Test_human_diff$Relationship)
  Test_human_diff_p$Relationship = Perm
  Test_human_diff_p %>% aov(distance ~ Relationship ,.) -> R2
  summary(R2)[[1]][4][[1]][1] -> F_value_p
  Perm_Fs = c(Perm_Fs, F_value_p)
}
sum(Perm_Fs >= F_value) / N_permutation
#Test difference between human-human and human-food
Distances_to_plot %>% filter(Distance %in% c("Human_to_Food", "Human_to_Human") ) %>% select(Distance, distance) -> Test_human_diff 
Test_human_diff %>% t.test(distance ~ Distance ,.) -> R_test2
R_test2$statistic -> F_value2
Perm_Fs2 = c()
set.seed(22311)
for (i in seq(1,N_permutation)){
  Test_human_diff_p = Test_human_diff 
  Perm = sample(Test_human_diff$Distance)
  Test_human_diff_p$Distance = Perm
  Test_human_diff_p %>% t.test(distance ~ Distance ,.) -> R2
  R2$statistic -> F_value_p
  Perm_Fs2 = c(Perm_Fs2, F_value_p)
}
sum(Perm_Fs2 <= F_value2) / N_permutation

#Test difference between mouse-mouse and human-human
Distances_to_plot %>% filter(Mice_1 == T & Mice_2 == T) %>% mutate(Distance = "Mouse-Mouse") %>% select(distance, Distance ) -> M_dist
Distances_to_plot %>% filter(Distance == "Human_to_Human") %>% select(distance, Distance ) -> H_dist
rbind(M_dist, H_dist) -> Test_human_diff
Test_human_diff %>% aov(distance ~ Distance ,.) -> R_test3
Test_human_diff %>% t.test(distance ~ Distance ,.) -> R_test3
R_test3$statistic -> F_value3
Perm_Fs3 = c()
set.seed(22311)
for (i in seq(1,N_permutation)){
  Test_human_diff_p = Test_human_diff 
  Perm = sample(Test_human_diff$Distance)
  Test_human_diff_p$Distance = Perm
  Test_human_diff_p %>% t.test(distance ~ Distance ,.) -> R2
  R2$statistic -> F_value_p
  Perm_Fs3 = c(Perm_Fs3, F_value_p)
}
sum(Perm_Fs3 <= F_value3) / N_permutation

#statistical test
Estimate_P = function(df, comparison = c("Same study", "Different continent") ){
  
  

}


Distances_to_plot$Relationship %>% unique()



Phenotypes_ba[ match( rownames(Distances) , Phenotypes_ba$ID_anal ), ] -> Phenotypes_ba
Phenotypes_ba %>% mutate(Group = ifelse(Reference_food == T, "Food", ifelse(Mice == T, "Mice", "Human"  )  ) ) %>% mutate(Group = ifelse( is.na(Group), "Mice", Group )) -> Phenotypes_ba
pairwiseAdonis::pairwise.adonis( as.dist(Distances) , Phenotypes_ba$Group , perm = 2000) -> Stats
#pairs Df    SumsOfSqs      F.Model           R2      p.value p.adjusted sig
#1 Human vs Food  1 3.610614e-03 1.639914e-03 3.031255e-06 0.7301349325 1.00000000    
#2 Human vs Mice  1 9.528016e+03 2.304985e+03 8.073545e-01 0.0004997501 0.00149925   *
#3  Food vs Mice  1 8.855826e+02 7.363584e+00 4.499982e-01 0.1799100450 0.53973013    
cophenetic.phylo( drop.tip(Bifido_animalis3, filter(Phenotypes_ba, Group != "Human" )$ID_anal ) ) -> Distances_human
adonis2( Distances_human ~ Continent + Continent:Country + Continent:Country:study_name , Phenotypes_ba[ match( rownames(Distances_human) , Phenotypes_ba$ID_anal ), ]  )
#Df SumOfSqs      R2      F Pr(>F)
#Continent                      2     0.61 0.00051 0.1244  0.400
#Continent:Country             17    23.03 0.01933 0.5547  0.283
#Continent:Country:study_name  30    92.86 0.07797 1.2676  0.175
#Residual                     440  1074.42 0.90218              
#Total                        489  1190.92 1.00000              


Distances_to_plot %>% mutate( Distance = ifelse(Distance == "Human_to_Human", Relationship, Distance ) ) %>% aov(distance ~ Distance, data = .) -> anova_results
TukeyHSD(anova_results)-> tukey_results
tukey_results$Distance %>% as.data.frame() %>% rownames_to_column("Comparison") %>% as_tibble() %>% arrange(`p adj`) %>% mutate(Sign = ifelse(`p adj`<0.05, T, F ) )
#Everything is sign different than Mice
#Nothing is sign diff than food
#Between countries, only same_study vs Different_continent is not sign. 
#these stats are inflated anyhow. Would need permutations if we need to support with stats

pairwise.adonis_s <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=999)
{
  
  co <- combn(unique(as.character(factors)),2)
  print(co)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  
  for(elem in 1:ncol(co)){
    print(elem)
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      }
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    x2 = data.frame(Fac = factors[factors %in% c(co[1,elem],co[2,elem])])
    
    ad <- adonis2(x1 ~ Fac, data = x2,
                  permutations = perm);
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$Df[1])
    SumsOfSqs <- c(SumsOfSqs,ad$SumOfSqs[1])
    F.Model <- c(F.Model,ad$F[1]);
    R2 <- c(R2,ad$R2[1]);
    p.value <- c(p.value,ad$`Pr(>F)`[1])
  }
  return(  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,)
)
}
#######################################
####FEATURES vs GEOGRAPHY EFFECT
######################################

read_tsv("/mnt/project/Make_Associations/Association/Results/Geography_cor/Association_results_merged.tsv") -> GeoEffect
#Genome length
read_tsv("/mnt/project/Make_Associations/Genome_characteristcs/Average_genome_length.tsv")  -> N_length
left_join(GeoEffect, N_length) %>% drop_na() %>% lm(Rho ~ avg_nt_length, .) %>% summary() -> summary_stats
write_tsv(summary_stats$coefficients %>% as.data.frame() %>% rownames_to_column("Feature") , "/mnt/project/Make_Associations/Genome_characteristcs/Results/GenomeLength_summary_stats.tsv")
#Traitar
read_tsv("/mnt/project/Make_Associations/Phenotypes/traitar_output_646SGBs_tab.txt") -> Traitar
read_tsv("/mnt/project/Make_Associations/Genome_characteristcs/traitar_output2.txt") %>% rename(SGB = sgb) -> Traitar2
missing_columns <- setdiff(names(Traitar), names(Traitar2))
for (col in missing_columns) { Traitar2[[col]] <- NA }
Traitar <- rbind(Traitar, Traitar2)

left_join(GeoEffect, Traitar) %>% drop_na() ->  Traitar_merged
Traitar_results = tibble()
for (Ann in colnames(Traitar) ){
  Traitar_merged -> Traitar_merged2
  Traitar_merged2[Traitar_merged2[Ann] == 0.5, Ann] = NA
  table( as.vector(Traitar_merged2[Ann]) )  -> Table
  if (dim(Table) == 1) { next }
  if (Table[1] < 20 || Table[2]<20 ){ next }
  if (Ann == "SGB"){ next }
  print(Ann)
  Formula = as.formula(paste0( "Rho ~ `", Ann, "`" ))
  lm(Formula, Traitar_merged2 %>% drop_na()) %>% summary() -> Model
  if (dim(Model$coefficients)[1] < 2) { next }
  Model$coefficients[ 2 , ] %>% t() %>% as_tibble() %>% mutate(Annotation = Ann) %>% rbind(Traitar_results, . ) -> Traitar_results
}
Traitar_results %>% arrange(`Pr(>|t|)`) %>% mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr") ) -> Traitar_results
Traitar_results %>%
write_tsv("/mnt/project/Make_Associations/Genome_characteristcs/Results/Traitar_summary_stats.tsv")

Traitar_merged %>% ggplot(aes(x= as.factor(`Casein hydrolysis`), y=Rho)) + geom_boxplot() + geom_sina() + theme_bw()


#include family
SGB_taxonomy = read_tsv("/mnt/project/Make_Associations/Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2", col_names = F ) %>% rename(SGB=X1, Taxonomy=X2)
SGB_taxonomy <- SGB_taxonomy %>% separate(Taxonomy, into = c("domain","phylum", "class", "order", "family", "genera", "species"), sep = "\\|", remove = FALSE)
SGB_taxonomy %>% mutate(SGB = paste0("t__", SGB  ) ) -> SGB_taxonomy
Traitar_merged %>% left_join(SGB_taxonomy %>% select(SGB, family) ) %>% mutate(Is_lachno = ifelse(family == "f__Lachnospiraceae", T, F) ) -> Traitar_merged
Traitar_results2 = tibble()
Traitar_results3 = tibble()
for (Ann in filter(Traitar_results, FDR<0.05)$Annotation ){
  Traitar_merged -> Traitar_merged2
  Traitar_merged2[Traitar_merged2[Ann] == 0.5, Ann] = NA
  if (Ann == "SGB"){ next }
  #Check if trait is present in both Lachno and not Lachno
  Traitar_merged2 %>% drop_na()  %>% group_by(!!sym(Ann),Is_lachno) %>% summarise(count = n()  )  -> Lachnoss
  if (dim(Lachnoss)[1] != 4){ cat("skipping", Ann) ;next }
  if ( sum(Lachnoss$count == 0) >0  ){ cat("skipping", Ann) ;next }
  print(Ann)
  Formula = as.formula(paste0( "Rho ~ Is_lachno +  `", Ann, "`" ))
  Formula2 =  as.formula(paste0( "Rho ~ Is_lachno *  `", Ann, "`" ))
  lm(Formula, Traitar_merged2 %>% drop_na()) %>% summary() -> Model
  lm(Formula2, Traitar_merged2 %>% drop_na()) %>% summary() -> Model2
  if (dim(Model$coefficients)[1] < 2) { next }
  Model$coefficients[ 3 , ] %>% t() %>% as_tibble() %>% mutate(Annotation = Ann) %>% rbind(Traitar_results2, . ) -> Traitar_results2
  Model2$coefficients[ 2:4 , ] %>% as.data.frame() %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Annotation = Ann) %>% rbind(Traitar_results3, . ) -> Traitar_results3
  
}
Traitar_results2 %>% mutate(FDR = p.adjust( `Pr(>|t|)`, "fdr" ) ) %>% write_tsv("/mnt/project/Make_Associations/Genome_characteristcs/Results/Traitar_summary_stats_lachno.tsv")
Traitar_results2 %>% select(Estimate, Annotation) %>% left_join(. , Traitar_results %>% select(Estimate, Annotation), by="Annotation", suffix=c("_with", "_without")  ) %>% mutate(Diff = Estimate_with - Estimate_without  )

library(ggrepel)
Traitar_results %>% ggplot(aes(x= Estimate , y= -log10(`Pr(>|t|)`), col=FDR<0.05 )) + geom_point()  + theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + 
   scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "#d62728")) + theme( axis.title = element_text(size = 16),  axis.text = element_text(size = 14)  ) + 
geom_text_repel( aes(label = ifelse(FDR < 0.05, as.character(Annotation), "")), box.padding = 0.5, point.padding = 0.5,  size = 5,   force = 5,  box.color = NA, direction = "both"  ) -> VolcanoTraitar
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/VolcanoTraitar_vs_GeographyRho.tiff",VolcanoTraitar)

library(ggforce)
To_plot =  c("Alkaline phosphatase", "Gram positive", "Trehalose", "Motile", "L-Rhamnose", "Catalase", "Casein hydrolysis", "Salicin", "Glucose fermenter")
Traitar_merged %>% select( c(Rho, FDR, SGB, To_plot) ) %>% gather(Trait, Available, c(To_plot) ) %>% filter(! Available == 0.5 ) %>% left_join(SGB_taxonomy %>% select(SGB, family) ) %>% mutate(Lachnospiraceae = ifelse(family == "f__Lachnospiraceae", T, F) ) %>%
  mutate(Available = ifelse(Available == 0, "No", "Yes")  ) %>% mutate(Trait = ifelse(Trait == "Alkaline phosphatase", "Alkaline\nphosphatase", ifelse(Trait == "Casein hydrolysis", "Casein\nhydrolysis", ifelse(Trait == "Glucose fermenter", "Glucose\nfermenter", Trait)) ) ) %>% 
  ggplot(aes(y=Rho, x= as.factor(Available), col=Lachnospiraceae  ))  + geom_sina(alpha=0.5) + geom_boxplot(outlier.shape = NA) + theme_bw() + facet_wrap(~Trait ) + scale_color_manual(values= c("FALSE" = "#370031", "TRUE" = "#CE8964" ) ) + xlab("Predicted trait") +
  theme(
    text = element_text(size = 11),  # Set the text size for all elements
    plot.title = element_text(size = 20, face = "bold"),  # Set title text size
    axis.title = element_text(size = 18, face = "bold"),  # Set axis title text size
    axis.text = element_text(size = 14),  # Set axis text size
    axis.ticks = element_line(size = 1.5),  # Set axis tick size
    legend.text = element_text(size = 14),  # Set legend text size
    legend.title = element_text(size = 16, face = "bold")  # Set legend title text size
  ) -> Trait_lachnospiraceae
library(patchwork)
VolcanoTraitar + Trait_lachnospiraceae -> Geography_ass_trait
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/GeographyRho_traits_composed.pdf", Geography_ass_trait, width = 11, height = 5.5)


Traitar_results2 %>% select(Estimate, Annotation) %>% left_join(. , Traitar_results %>% select(Estimate, Annotation), by="Annotation", suffix=c("_with", "_without")  ) %>% mutate(Percentage_change = 1 - (Estimate_with/Estimate_without  ) ) %>% summarise(mean(Percentage_change))

#For the traits that are exclusive of Lachno o no exclusive of lachno
Traitar_results4 = tibble()
for (Ann in filter(Traitar_results, FDR<0.05)$Annotation ){
  if (Ann %in% Traitar_results2$Annotation){ next }
  Traitar_merged -> Traitar_merged2
  Traitar_merged2[Traitar_merged2[Ann] == 0.5, Ann] = NA
  Traitar_merged2 %>% drop_na()  %>% group_by(!!sym(Ann),Is_lachno) %>% summarise(count = n()  )  -> Lachnoss
  if ( dim( filter(Lachnoss, Is_lachno == T) )[1] < 2  ){ With_lachno = F ; Traitar_merged2 %>% filter(Is_lachno == F ) -> Traitar_merged2  
  } else { With_lachno = T; Traitar_merged2 %>% filter(Is_lachno == T ) -> Traitar_merged2  }
  Formula = as.formula(paste0( "Rho ~  `", Ann, "`" ))
  lm(Formula, Traitar_merged2 %>% drop_na()) %>% summary() -> Model
  Model$coefficients[ 2 , ] %>% t() %>% as_tibble() %>% mutate(Annotation = Ann) %>% mutate(With_Lachnospiraceae = With_lachno) %>% rbind(Traitar_results4, . ) -> Traitar_results4
}  
Traitar_results4 %>% select(Estimate, Annotation) %>% left_join(. , Traitar_results %>% select(Estimate, Annotation), by="Annotation", suffix=c("_with", "_without")  ) %>% mutate(Diff = Estimate_with - Estimate_without  ) %>%mutate(Percentage_change = 1 - (Estimate_with/Estimate_without  ) ) %>% summarise(mean(abs(Percentage_change)))

  

#Prevalence
read_tsv("/mnt/project/Make_Associations/Genome_characteristcs/prevalence_SGB_environments.tsv") -> Prev_env
left_join(GeoEffect, Prev_env) %>% drop_na() -> Prev_env
Prev_env %>% group_by(SGB) %>% summarise(Prevalence_total = sum(prevalence_bin), Rho )  %>% distinct(SGB, .keep_all = T) %>% lm(Rho ~ Prevalence_total, .) %>% summary()
Environments_results = tibble()
Environments_results_proportion = tibble()

for (Env in Prev_env$environment %>% unique()){
  Prev_env %>% filter( environment==Env ) -> ForModel
  lm(Rho ~ prevalence_bin, ForModel ) %>% summary() -> Model
  lm(Rho ~ rel_prevalence, ForModel ) %>% summary() -> Model2
  
  if (dim(Model$coefficients)[1] < 2) { next }
  Model$coefficients[ 2 , ] %>% t() %>% as_tibble() %>% mutate(Environment = Env) %>% rbind(Environments_results, . ) -> Environments_results
  Model2$coefficients[ 2 , ] %>% t() %>% as_tibble() %>% mutate(Environment = Env) %>% rbind(Environments_results_proportion, . ) -> Environments_results_proportion
  
}
Environments_results %>% arrange(`Pr(>|t|)`) %>% mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr") ) %>% mutate(Model = "Binary") -> Environments_results
Environments_results_proportion %>% arrange(`Pr(>|t|)`) %>% mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr") ) %>% mutate(Model = "Proportion") -> Environments_results_proportion
rbind(Environments_results, Environments_results_proportion) %>% write_tsv("/mnt/project/Make_Associations/Genome_characteristcs/Results/EnvironmentPrevalence.tsv")

Prev_env %>% filter(environment %in% c("Stool_Ancient", "Wild_Mice") ) %>%     ggplot(aes(x= as.factor(prevalence_bin), y = Rho )) + geom_boxplot(outlier.shape = NA) + theme_bw() + facet_wrap(~environment) + geom_sina()
Prev_env %>% filter(environment %in% c("Stool_Ancient", "Stool_NW", "NHP") ) %>%     ggplot(aes(x= rel_prevalence, y = Rho )) + geom_point() + theme_bw() + facet_wrap(~environment) + geom_smooth(method="lm")
#stool not westernized (Stool_NW) , NHP (non-human primate)


#Heatmao
pivot_wider(Prev_env %>% select(SGB, environment, rel_prevalence) , names_from = environment, values_from = rel_prevalence) %>% left_join(. ,Prev_env %>% select(SGB, Rho)  ) %>% select(-SGB) -> Matrix_pheatmap
data.frame(Label = colnames(Matrix_pheatmap) ) -> col_labels
col_colors <- c(rep("black", 16), "red")
names(col_colors) = col_labels$Label
data.frame( Feature = c(rep("Prevalence", 16), "Correlation coefficient") ) -> Ann
rownames(Ann) = colnames(Matrix_pheatmap)
Matrix_pheatmap %>% cor() %>% pheatmap::pheatmap( annotation_col =   Ann  )


pivot_wider(Prev_env %>% select(SGB, environment, rel_prevalence, Rho) , names_from = environment, values_from = rel_prevalence)  %>% lm(Rho ~ Stool_Ancient + Stool_NW+ NHP , . ) %>% summary()


Prev_env %>% filter(environment %in% c("Stool_Ancient") ) %>% mutate(prevalence_bin = ifelse(prevalence_bin=="0", "SGB not present", "SGB present") ) %>%  ggplot(aes(x= as.factor(prevalence_bin) , y = Rho )) + 
  geom_boxplot(outlier.shape = NA) + geom_sina(alpha=0.5) + theme_bw() + ylab("Geographical effect (rho)") + 
  xlab("SGB present in Ancient stool samples") +theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) -> AncientBoxplot
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/SGBinAncient_vs_GeogrpahyRho.tiff",AncientBoxplot, width = 5.6, height = 5.6)






###Sup fig 2: show distributions related vs unrelated
Prefix = "/mnt/project/"

read_csv(paste0(Prefix, "Make_Associations/Association/Results/Summaries/Protocol/SequencingProtocol_results.csv") ) -> Tested_SGB
read_tsv(paste0(Prefix, "Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv")) -> Metadata

Clean_names =  function(Name){
  'Samples from 500FG FSK need everything after _ to be removed. If _metaphlan4 is attached at the end of the name, remove'
  if ( grepl("HV", Name)) {
    Name = str_split(Name, "_")[[1]][1]
  }   
  Name = str_replace(Name, "_metaphlan4", "") 
  return(Name)
}

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

#Get tree
Get_distance_dist = function(SGB){
  Tree = paste0(paste0(Prefix,"Symlink_phylo/IQtree.",SGB, ".TreeShrink.tre"))
  if (! file.exists(Tree)){ Tree =paste0(paste0(Prefix,"Symlink_phylo/RAxML_bestTree.",SGB, ".TreeShrink.tre"))  }
  Process_tree(Tree) -> Tree
  Metadata %>% filter(study_name %in% c("500FG_FSK", "SchirmerM_2016") ) -> Metadata_rep
  Metadata_rep %>% filter(ID_anal %in% Tree$tip.label ) %>% select(study_name, subject_id, ID_anal) -> Metadata_rep
  Metadata_rep %>% group_by(subject_id) %>% summarise(N = n()) %>% filter(N > 1) -> Keep
  Metadata_rep %>% filter(subject_id %in% Keep$subject_id) -> Metadata_rep
  keep.tip(Tree, Metadata_rep$ID_anal) -> Tree
  #Get distribution of distances
  cophenetic.phylo(Tree) -> Distances
  dist_df <- as.data.frame(as.table(Distances))
  colnames(dist_df) <- c("rowname", "colname", "distance")
  as_tibble(dist_df) %>% filter(!rowname==colname) -> dist_df
  #Add info
  left_join(dist_df, Metadata_rep, by=c("rowname" = "ID_anal") ) %>% left_join(., Metadata_rep, by=c("colname" = "ID_anal"),  suffix=c("_1", "_2" ) ) -> dist_df
  dist_df %>% filter(! rowname == colname   ) -> dist_df
  dist_df %>% mutate(Same_subject =  ifelse( subject_id_1 == subject_id_2, T, F ) ) -> dist_df
  sd(dist_df$distance) -> Norm_factor
  dist_df %>% mutate(distance = distance/Norm_factor) -> dist_df
  #dist_df %>% ggplot(aes(x=distance, fill=Same_subject )) + geom_histogram() + theme_bw() %>% print()
  return(dist_df)
}
All_distances = tibble()
for (SGB in unique(Tested_SGB$SGB) ){
  Get_distance_dist(SGB) %>% mutate(SGB= SGB) %>% rbind(All_distances, . )  -> All_distances
}


All_distances %>% ggplot(aes(x=distance, fill=Same_subject )) + geom_histogram() + theme_bw() 
All_distances %>% ggplot(aes(y=distance, x=Same_subject )) + geom_boxplot(outlier.size = 0.5 ) + theme_bw() + coord_flip() + ylab("Phylogenetic distance") + xlab("Sample pair from same individual") + 
  theme(
    text = element_text(size = 16),  # Set the text size for all elements
    plot.title = element_text(size = 20, face = "bold"),  # Set title text size
    axis.title = element_text(size = 18, face = "bold"),  # Set axis title text size
    axis.text = element_text(size = 14),  # Set axis text size
    axis.ticks = element_line(size = 1.5),  # Set axis tick size
    legend.text = element_text(size = 14),  # Set legend text size
    legend.title = element_text(size = 16, face = "bold")  # Set legend title text size
  )


All_distances %>% filter(Same_subject == T) %>% arrange(desc(distance)) %>% filter(distance > 1) %>% group_by(SGB) %>% summarise(n())

#test with permutations
All_distances %>% select(distance, Same_subject) -> To_Test
To_Test %>% t.test(distance ~ Same_subject, . ) -> Original
Original$statistic -> orginal_t
Permutation_n = 1000
Permuted_t = c()
set.seed(1341341)
for (i in seq(1, Permutation_n)){
  To_Test_p = To_Test
  sample(To_Test$Same_subject) -> To_Test_p$Same_subject
  To_Test_p %>% t.test(distance ~ Same_subject, . ) -> Permuted
  Permuted_t = c(Permuted_t, Permuted$statistic)
}
sum(Permuted_t >= orginal_t)/ Permutation_n

#Test: Avg FSK vs APK same sample / avg( (FSG vs FSK) (APK vs APK) )
#Get info abundance
read_rds("/mnt/project/SGB_abundances/Data/Merged_Groningen_SGB.rds") -> Abundance
colnames(Abundance) %>% sapply(function(x){ if(grepl("\\|", x) ){ str_split(x, "\\|") -> y ; y[[1]] -> y ;  y[[length(y)]] } else{ x }  } ) -> NewNames
colnames(Abundance) = NewNames
str_replace(Abundance$ID, "_APK", "") -> Abundance$ID
Abundance %>% filter(ID %in% unique(All_distances$rowname) ) -> Abundance

Ratios_tibble = tibble()
for (S in unique(All_distances$SGB) ){
  #Get distances genetics in tree
  All_distances %>% filter(SGB == S) -> Distances_SGB
  #Get distances abundance in tree
  Abundance %>% select("ID", S) -> Abundance_SGB
  Abundance_SGB %>% as.data.frame() %>% column_to_rownames("ID") %>% dist("euclidean") ->Dist_abundance
  dist_df <- as.data.frame(as.table(as.matrix(Dist_abundance)))
  colnames(dist_df) <- c("rowname", "colname", "distance_ab")
  as_tibble(dist_df) %>% filter(!rowname==colname) -> dist_df
  Distances_SGB %>% left_join(dist_df) -> Distances_SGB
  #Prepare for ratios
  Distances_SGB %>% mutate(Same_study = ifelse(study_name_1 == study_name_2, T, F) ) -> Distances_SGB
  Distances_same = Distances_SGB %>% filter(Same_subject == T)
  Distance_diff_dstudy = Distances_SGB %>% filter(Same_subject == F & Same_study == F)
  #Ratios genetics
  Numerator = mean(Distances_same$distance)
  Denominator = mean(Distance_diff_dstudy$distance)
  Ratio_tree = Numerator / Denominator
  #Ratios abundance
  Numerator = mean(Distances_same$distance_ab)
  Denominator = mean(Distance_diff_dstudy$distance_ab)
  Ratio_ab = Numerator / Denominator
  Ratios_tibble = rbind(Ratios_tibble, tibble(SGB = S, Ratio_genetics=Ratio_tree, Ratio_Abundance=Ratio_ab) )
}  
Ratios_tibble %>% drop_na() -> Ratios_tibble
Ratios_tibble %>% gather(Ratio, Distance, c("Ratio_genetics", "Ratio_Abundance") ) %>% mutate(Ratio = ifelse(Ratio == "Ratio_genetics", "Genetic", "Abundance" ) ) ->Ratios_tibble2
Ratios_tibble2  %>% ggplot(aes(x=Ratio, y=Distance))  + geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.6) + 
  geom_line(aes(group = SGB), alpha = 0.2, linewidth=0.2) + theme_bw() + xlab("Ratio average same subject distance \n average distance between unrelated APK-FSK samples") + coord_flip() -> Plot_ratioDistance
wilcox.test( Ratios_tibble$Ratio_Abundance, Ratios_tibble$Ratio_genetics, paired = T )
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/FSK_vs_APK_distanceratios.pdf",Plot_ratioDistance, width = 4, height=4 )

#rectale tre
Process_tree("/mnt/project/Symlink_phylo/IQtree.t__SGB4933_group.TreeShrink.tre") -> rectale
Metadata %>% filter(ID_anal %in% rectale$tip.label ) -> Metadata_rectale
Metadata_rectale %>% distinct(ID_anal, .keep_all = T ) -> Metadata_rectale

set.seed(823)
sample(Metadata_rectale$ID_anal, 1000) -> random_subsample
Metadata_rectale %>% filter(ID_anal %in% random_subsample ) -> Metadata_rectale2
keep.tip(rectale, random_subsample) -> rectale2
 
Metadata_rectale2 %>% mutate(Study = ifelse(study_name == "DAG3", "GacesaR_2023", ifelse(study_name == "LLD", "ZhernakovaA_2016", ifelse(study_name %in% c("HMP_2019_ibdmdb", "AsnicarF_2021", "MehtaRS_2018", "ZeeviD_2015", "YachidaS_2019" ), study_name, "Other" )  )  )  ) -> Metadata_rectale2

ggtree(rectale2, layout = "fan",  open.angle=15, size=0.1 ) %<+% select(Metadata_rectale2, c(ID_anal, Continent, Country, Study ) ) +
geom_tiplab(aes(col=Continent, label=Country), size=1 ) + scale_color_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) +
new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) +
new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Study), width=0.03,offset=0.1 ) + scale_fill_manual(values=c("GacesaR_2023" = "#6A3D9A", "ZhernakovaA_2016" = "#CAB2D6", "HMP_2019_ibdmdb" = "#FB9A99", "AsnicarF_2021" ="skyblue2", "MehtaRS_2018" = "yellow4", "ZeeviD_2015" = "brown", "YachidaS_2019" = "darkorange4", Other="grey"  ))


c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")



##MAke plot for bacteria with multiple assocations / Blautia_wexlerae

set.seed(823)
Europe_filtered %>% filter(SGB == "t__SGB4837_group")
Metadata %>% filter(Continent == "Europe") %>%  select(ID_anal, study_name, Age, IBD, Mesalazines, CalciumSupplements, VitB12, melanoma, PPI, VitD) -> Meta_for_tree
Process_tree("/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB4837_group.TreeShrink.tre") -> Tree
#Meta_for_tree %>% filter(ID_anal %in% Tree$tip.label )  %>% distinct(ID_anal, .keep_all = T) -> Meta_for_tree
#keep.tip(Tree, Meta_for_tree$ID_anal) -> Tree

Meta_for_tree = Meta_for_tree %>% mutate(Age_group = ifelse(Age ==0, "<1", ifelse(Age < 15, "(0-15)", ifelse( Age < 30, "[15, 30)", ifelse( Age < 70, "[30, 70)", ">=70") )))) %>% mutate(Age_group = factor(Age_group, levels=c( "<1", "(0-15)", "[15, 30)", "[30, 70)", ">=70"  )))
Meta_for_tree %>% filter(! is.na(Age) ) -> Meta_for_tree
Meta_for_tree %>% filter(ID_anal %in% Tree$tip.label )  %>% distinct(ID_anal, .keep_all = T) -> Meta_for_tree
keep.tip(Tree, Meta_for_tree$ID_anal) -> Tree

sample(Meta_for_tree$ID_anal, 500) -> random_subsample


filter(Meta_for_tree,IBD == 1)$ID_anal %>% c(random_subsample, . ) -> random_subsample
filter(Meta_for_tree, melanoma == 1)$ID_anal %>% c(random_subsample, . ) -> random_subsample
filter(Meta_for_tree, PPI == 1)$ID_anal %>% c(random_subsample, . ) -> random_subsample
unique(random_subsample) -> random_subsample


Meta_for_tree %>% filter(ID_anal %in% random_subsample ) -> Meta_for_tree2
keep.tip(Tree, random_subsample) -> Tree2



ggtree(Tree2, layout = "fan",  open.angle=15, size=0.1 ) %<+% Meta_for_tree2 +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age_group), width=0.03,offset=0  ) + scale_fill_manual(values = c("<1" = "#5E1DB5",  "(0-15)"= "#5A3C82", "[15, 30)" = "#0E27E8", "[30, 70)"= "#EBB649",  ">=70"="#B5731D" ), na.value="white" ) +   labs(fill = "Age") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(IBD)), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "#d62728", "0" = "grey"), na.value ="white", labels = c("No", "Yes") ) + labs(fill = "IBD") +
#  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(Mesalazines) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "#1f77b4", "0" = "grey"), na.value ="white" ) + 
#  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(CalciumSupplements) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "gold1", "0" = "grey"), na.value ="white" ) +
#  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(VitB12) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "palegreen2", "0" = "grey"), na.value ="white" ) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(melanoma) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "black", "0" = "grey"), na.value ="white", labels = c("No", "Yes") ) + labs(fill = "Melanoma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(PPI) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "brown", "0" = "grey"), na.value ="white", labels = c("No", "Yes") ) + labs(fill = "PPI") -> p
#  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(VitD) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "#FF7F00", "0" = "grey"), na.value ="white" ) 
  
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/Blautia_multipleAssociations.pdf", p)  
  


#MAke plot for bacteria with multiple assocations: Collinsella
Europe_filtered %>% filter(SGB == "t__SGB14546_group") 
Metadata %>% filter(Continent == "Europe") %>%  select(ID_anal, study_name, Age, IBD, Hypertension, Asthma, Sex, melanoma, CrohnsDisease, UlcerativeColitis) -> Meta_for_tree
Process_tree("/mnt/project/Symlink_phylo/IQtree.t__SGB14546_group.TreeShrink.tre") -> Tree
Meta_for_tree = Meta_for_tree %>% mutate(Age_group = ifelse(Age ==0, "<1", ifelse(Age < 15, "(0-15)", ifelse( Age < 30, "[15, 30)", ifelse( Age < 70, "[30, 70)", ">=70") )))) %>% mutate(Age_group = factor(Age_group, levels=c( "<1", "(0-15)", "[15, 30)", "[30, 70)", ">=70"  )))
Meta_for_tree %>% filter(! is.na(Age) ) -> Meta_for_tree
Meta_for_tree %>% filter(ID_anal %in% Tree$tip.label )  %>% distinct(ID_anal, .keep_all = T) -> Meta_for_tree
keep.tip(Tree, Meta_for_tree$ID_anal) -> Tree


set.seed(13413)
sample(Meta_for_tree$ID_anal, 500) -> random_subsample
filter(Meta_for_tree,IBD == 1)$ID_anal %>% c(random_subsample, . ) -> random_subsample
filter(Meta_for_tree, melanoma == 1)$ID_anal %>% c(random_subsample, . ) -> random_subsample
filter(Meta_for_tree, Hypertension == 1)$ID_anal %>% c(random_subsample, . ) -> random_subsample
filter(Meta_for_tree, Asthma == 1)$ID_anal %>% c(random_subsample, . ) -> random_subsample
unique(random_subsample) -> random_subsample

Meta_for_tree %>% filter(ID_anal %in% random_subsample ) -> Meta_for_tree2
keep.tip(Tree, random_subsample) -> Tree2

ggtree(Tree2, layout = "fan",  open.angle=15, size=0.1 ) %<+% Meta_for_tree2 +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age_group), width=0.03,offset=0  ) + scale_fill_manual(values = c("<1" = "#5E1DB5",  "(0-15)"= "#5A3C82", "[15, 30)" = "#0E27E8", "[30, 70)"= "#EBB649",  ">=70"="#B5731D" ), na.value="white" ) +   labs(fill = "Age") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(IBD)), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "#d62728", "0" = "grey"), na.value ="white", labels = c("No", "Yes") ) + labs(fill = "IBD") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(Sex) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "#1f77b4", "0" = "grey"), na.value ="white" ) + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(Hypertension) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "gold1", "0" = "grey"), na.value ="white" ) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(Asthma) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "lightblue", "0" = "grey"), na.value ="white" ) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(melanoma) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "black", "0" = "grey"), na.value ="white", labels = c("No", "Yes") ) + labs(fill = "Melanoma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(UlcerativeColitis) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "brown", "0" = "grey"), na.value ="white", labels = c("No", "Yes") ) + labs(fill = "UC") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=as.factor(CrohnsDisease) ), width=0.03,offset=0.1  ) + scale_fill_manual(values=c("1" = "darkred", "0" = "grey"), na.value ="white", labels = c("No", "Yes") ) + labs(fill = "CD") -> p



###############################
###B. longum and infants#######
###############################
set.seed(1422)


Remove_repeated = function(Metadata, pruned_tree){
  'Check which samples are repeated (either duplicates or longitudinally) and only picks one from the repeated'
  set.seed(89777)
  Metadata %>% group_by(subject_id) %>% summarise(N = n()) %>% arrange(desc(N)) %>% filter(N > 1 ) -> Replicates
  remove = tibble()
  for (i in Replicates$subject_id){
    Metadata %>% filter(subject_id == i) -> Filtered ; Samples = Filtered$ID_anal
    Keep = sample(Samples, 1)
    remove = c(remove, Samples[!Samples==Keep] ) 
  }
  remove = unlist(remove)
  Metadata %>% filter(! ID_anal %in% remove) -> Metadata
  pruned_tree %>% drop.tip(remove) -> pruned_tree
  return(list(Metadata, pruned_tree))
}
sort_values <- function(value1, value2) {
  if (value1 < value2) {
    return(paste0(value1, value2))
  } else {
    return(paste0(value2, value1))
  }
}

Tree = Process_tree("/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB17248.TreeShrink.tre")
Metadata %>% filter(ID_anal %in% Tree$tip.label ) %>% Remove_repeated(., Tree) -> Out
Out[[1]] -> Meta_bifido
Out[[2]] -> Tree_bifido

filter(Meta_bifido, infant == T)$ID_anal -> Keep_infant
filter(Meta_bifido %>% filter(Continent == "Africa"), infant == F)$ID_anal %>% sample(. , size=40 ) -> Keep_Adult
filter(Meta_bifido %>% filter(Continent == "Asia"), infant == F)$ID_anal %>% sample(. , size=40 ) %>% c(Keep_Adult, . ) -> Keep_Adult
filter(Meta_bifido %>% filter(Continent == "Europe"), infant == F)$ID_anal %>% sample(. , size=100 ) %>% c(Keep_Adult, . ) -> Keep_Adult
filter(Meta_bifido %>% filter(Continent == "North_America"), infant == F)$ID_anal %>% sample(. , size=40 ) %>% c(Keep_Adult, . ) -> Keep_Adult
filter(Meta_bifido %>% filter(Continent == "South_America"), infant == F)$ID_anal %>% sample(. , size=40 ) %>% c(Keep_Adult, . ) -> Keep_Adult

Meta_bifido %>% filter( ID_anal %in%  c(Keep_Adult, Keep_infant) ) -> Meta_bifido2 
keep.tip(phy = Tree_bifido, tip=Meta_bifido2$ID_anal) -> Tree_bifido2

Distances = cophenetic.phylo(Tree_bifido2) 
dist_df <- as.data.frame(as.table(Distances))
colnames(dist_df) <- c("rowname", "colname", "distance")
as_tibble(dist_df) %>% filter(!rowname==colname) -> dist_df
#Add info
left_join(dist_df, Meta_bifido2 %>% select(ID_anal, Country, Continent, infant ), by=c("rowname" = "ID_anal") ) %>% left_join(., Meta_bifido2 %>% select(ID_anal, Country, Continent, infant ), by=c("colname" = "ID_anal"),  suffix=c("_1", "_2" ) ) -> dist_df
dist_df %>% filter(! rowname == colname   ) -> dist_df
sd(dist_df$distance) -> Norm_factor
dist_df %>% mutate(distance = distance/Norm_factor) -> dist_df
dist_df %>% mutate(Comparison = paste(pmin(infant_1, infant_2), pmax(infant_1, infant_2), sep="_")) -> dist_df
dist_df %>% mutate(Continent_Comparison = paste(pmin(Continent_1, Continent_2), pmax(Continent_1, Continent_2), sep="_")) -> dist_df
dist_df %>% filter(! Continent_1 %in% c("North_America", "South_America", "Asia" ) | Continent_2 %in% c("North_America", "South_America", "Asia" )  ) -> dist_df

dist_df  %>%
  ggplot(aes(x=Comparison, y=distance, col=Continent_Comparison )) + geom_boxplot(outlier.shape=NA) + geom_sina() + theme_bw()



######################################
####Check assocaiton collinsella######
######################################

Metadata %>% filter(study_name %in% c("FrankelAE_2017","GopalakrishnanV_2018","LeeKA_2022","MatsonV_2018","McCullochJA_2022", "PetersBA_2019", "WindTT_2020") ) -> meta_melanoma
#add age and other phenos
Tree_col = Process_tree("/mnt/project/Symlink_phylo/IQtree.t__SGB14546_group.TreeShrink.tre")
meta_melanoma %>% select(ID_anal, study_name, Country, Continent,  Age,melanoma ,disease) %>% filter(ID_anal %in% Tree_col$tip.label) %>% mutate(melanoma = melanoma) -> meta_melanoma
set.seed(123)
c(meta_melanoma$ID_anal,  sample(Tree_col$tip.label, 1000)) -> Subselection
keep.tip(Tree_col, Subselection) ->  Tree_col2

ggtree(Tree_col2, layout="fan", open.angle=15, size=0.1) %<+% meta_melanoma  -> p
p + geom_fruit( geom="geom_tile", mapping = aes(fill=melanoma), width=0.03,offset=0.1)  +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Country), width=0.03,offset=0.1 )  + #scale_fill_manual(values=c25)
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=study_name), width=0.03,offset=0.1 ) + scale_fill_manual(values=c25) +
  geom_balance(node=1364  , fill="skyblue", color=NA, alpha=0.3) #+ geom_balance(node=1366  , fill="red", color=NA, alpha=0.3)

#highlight clade
Names = filter(meta_melanoma, study_name == "WindTT_2020")$ID_anal
match( Names , Tree_col2$tip.label) -> CHECK
common_ancestor <- getMRCA(Tree_col2, CHECK)
Names = sample(filter(meta_melanoma, study_name == "LeeKA_2022")$ID_anal, 5)
print(Names)
match( Names , Tree_col2$tip.label) -> CHECK2
common_ancestor <- getMRCA(Tree_col2, CHECK2)
common_ancestor

Tree_col2$tip.label[getDescendants(Tree_col2, node = 1364)] -> Melanoma_clade
Melanoma_clade 
Metadata %>% filter(ID_anal %in% Subselection) %>% mutate(Meta_clade = ifelse(ID_anal %in% Melanoma_clade, 1, 0) ) -> Metadata_mel
Metadata_mel %>% mutate(melanoma = ifelse(is.na(melanoma), 0, melanoma )) %>% select(melanoma, Age, Meta_clade) %>% drop_na() %>% glm( melanoma ~ Age + Meta_clade,., family="binomial" ) %>% summary()

#Association with previous treatment

#read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/McCullochJA.tsv") %>% select(sample_id, subject_id, previous_therapy, treatment) -> MC_meta
#one_hot_encode <- function(df, column_name) {
  df %>%
    separate_rows({{column_name}}, sep = ";") %>%
    pivot_wider(names_from = {{column_name}}, values_from = {{column_name}}, values_fn = length, values_fill = 0)
}
#MC_meta %>% one_hot_encode(treatment) %>% rename(ID_anal = sample_id) -> MC_meta2
#MC_meta2 %>% left_join(., Metadata_mel, by="ID_anal")  %>% mutate(melanoma = ifelse(is.na(melanoma), 0, melanoma )) %>% select( one_of(c("melanoma", "Age", "Meta_clade", colnames(MC_meta2)[3:12] ))) %>% drop_na() %>%
#  glm( Meta_clade ~ . ,., family="binomial" ) %>% summary()

#read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/Leeka2022.tsv") %>% select(sample_id, subject_id, history_of_therapy) -> LK_meta
#LK_meta %>% rename(ID_anal = sample_id)  %>% left_join(., Metadata_mel, by="ID_anal")  %>% mutate(melanoma = ifelse(is.na(melanoma), 0, melanoma )) %>% select( one_of(c("melanoma", "Age", "Meta_clade", "history_of_therapy" ))) %>% drop_na() %>%
#  glm( Meta_clade ~ history_of_therapy ,., family="binomial" ) %>% summary()


read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/McCullochJA.tsv") %>% select(sample_id, subject_id, previous_therapy, treatment) %>% mutate(Clade = ifelse(sample_id %in% Metadata_mel$sample_id, T, F )) -> MC_meta
MC_meta %>% glm( Clade ~ previous_therapy ,., family="binomial" ) %>% summary()
read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/Leeka2022.tsv") %>% select(sample_id, subject_id, history_of_therapy) %>% mutate(Clade = ifelse(sample_id %in% Metadata_mel$sample_id, T, F )) ->LK_meta
LK_meta %>% glm( Clade ~ history_of_therapy ,., family="binomial" ) %>% summary()

MC_meta %>% select(Clade, previous_therapy) %>% mutate(Cohort = "MC") %>%
  rbind(LK_meta %>% mutate(previous_therapy=history_of_therapy, Cohort="LK") %>% select(Clade, previous_therapy, Cohort)  ) %>%
  glmer( Clade ~ previous_therapy + (1|Cohort)  ,., family="binomial" ) %>% summary()




###########################################
##Functional analysis in melanoma#####
###########################################
library(tidyverse)
library(anpan)
library(lmerTest)
c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")

Process_RDS = function(Data){
  read_rds(Data) -> Data
  Data$model_input -> Info
  Data$pglmm_fit$summary() -> Data
  print("Preparing information table")
  Data %>% filter( grepl("phylo", variable))  %>% filter(! grepl("std", variable)) %>% mutate(sample_id = factor(c("-", as.character(Info$sample_id)))) %>%
    select(sample_id, median) %>% rename(phylo_effect_median = median) -> Phylo
  left_join(Info, Phylo) %>% mutate(sample_id = as.character(sample_id)) -> Phylo
  return(Phylo)
}
one_hot_encode <- function(df, column_name) {
  df %>%
    separate_rows({{column_name}}, sep = ";") %>%
    pivot_wider(names_from = {{column_name}}, values_from = {{column_name}}, values_fn = length, values_fill = 0)
}
Process_tree = function(Tree){
  read.tree(Tree) -> Tree
  if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }
  Tree = midpoint.root(Tree)
  #2. Rename labels tree
  print("Processing tree")
  Tree$tip.label %>% sapply(. , Clean_names) -> New_names
  Tree$tip.label = as.vector(New_names)
  #if ( length(Tree$tip.label) < 100 ) { print("Tree is too small") ; q() }
  return(Tree)
  
}


Process_RDS("/mnt/project/Make_Associations/Association/Results/Continent_stratified/Age,Country_Europe/Models/t__SGB14546_group/melanoma/Model.rds") -> melanoma_europe
Process_RDS("/mnt/project/Make_Associations/Association/Results/Continent_stratified/Age,Country_North_America/Models/t__SGB14546_group/melanoma/Model.rds") -> melanoma_NA
#Merge North American and European results
rbind(melanoma_europe %>% mutate(Continent = "Europe") , melanoma_NA %>% mutate(Continent = "North_America") ) -> Phylo
#Add other covariates: Include desired covariates in the select statment
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> Phenos
Phenos %>% select(study_name, ID_anal, PC, Sex) %>% rename(sample_id= ID_anal) %>% left_join(Phylo, . ) -> Info
#Add study-specific covariates
read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/McCullochJA.tsv") %>% select(sample_id, subject_id, previous_therapy, treatment, RECIST) %>% rename(ID_anal = sample_id)  -> MC_meta
MC_meta %>% one_hot_encode(treatment)-> MC_meta2
read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/Leeka2022.tsv") %>% select(sample_id, subject_id, history_of_therapy, RECIST) %>% rename(ID_anal = sample_id)  -> LK_meta
read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/Franke2017.tsv") %>% select(sample_id, subject_id, RECIST) %>% rename(ID_anal = sample_id)  -> Franke_meta
Franke_meta$ID_anal %>% sapply(function(x){str_split(x, "_")[[1]][2] }) %>% as.vector() -> Franke_meta$ID_anal

Info %>% rename(ID_anal = sample_id) %>% left_join(MC_meta2, by="ID_anal") %>% left_join(LK_meta, by="ID_anal") %>% left_join(Franke_meta, by="ID_anal") -> Info
Info  %>% mutate(Study = ifelse(melanoma == 1|PC==1 , study_name, "control_study") )  -> Info
Info %>% mutate(RECIST = ifelse(is.na(RECIST.x) & is.na(RECIST.y) & is.na(RECIST), NA, ifelse(!is.na(RECIST.x), RECIST.x, ifelse( !is.na(RECIST.y), RECIST.y, RECIST  ) )) ) -> Info
Info %>% mutate(PC = ifelse(PC == 1, T, F) ) -> Info


#Read Tree
Tree = Process_tree("/mnt/project/Symlink_phylo/IQtree.t__SGB14546_group.TreeShrink.tre")
Phenos %>% select(ID_anal, study_name, melanoma, PC, Sex) %>% rename(sample_id= ID_anal) %>% filter(sample_id %in% Tree$tip.label ) %>% mutate(melanoma = ifelse(is.na(melanoma), 0 , melanoma ), PC=ifelse(is.na(PC), 0, PC ) )  -> InfoAll
keep.tip(Tree, InfoAll$sample_id ) -> Tree_raw
keep.tip(Tree, Info$ID_anal) -> Tree


#Select clade
Names = filter(Info, study_name == "WindTT_2020")$ID_anal
match( Names , Tree$tip.label) -> CHECK
common_ancestor <- getMRCA(Tree, CHECK)
getParent(Tree, 3763)
#clade members
Tree$tip.label[getDescendants(Tree, 3760)] -> Members_clade
Members_clade[!is.na(Members_clade)] -> Members_clade
Info %>% mutate(Melanoma_clade = ifelse(ID_anal %in% Members_clade, T, F) ) -> Info
Info %>% mutate(Country_sample = ifelse(Country %in% c("NLD", "GBR", "USA", "CHE"), Country, "Other" ) ) -> Info
#clade members overall
getMRCA(Tree_raw, match( Members_clade, Tree_raw$tip.label  ) ) #9608
Tree_raw$tip.label[getDescendants(Tree_raw, 9608)] -> Members_clade2
InfoAll %>% mutate(Melanoma_clade = ifelse(sample_id %in% Members_clade2, T, F) ) -> InfoAll
InfoAll %>% filter(Melanoma_clade == 1) %>% filter(study_name %in% c("LLD", "IBD", "DAG3") )

#Make plot
ggtree(Tree, layout="fan", open.angle=15, size=0.1) %<+% Info  -> p
p + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=melanoma), width=0.03,offset=0.1) + scale_fill_manual(values=c("TRUE"= "#E31A1C", "FALSE"="white") )  + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=PC), width=0.03,offset=0.1) + scale_fill_manual(values=c("TRUE"= "blue1", "FALSE"="white") )  + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Study), width=0.03,offset=0.1 ) +     scale_fill_manual(values=c("control_study"="grey", "FrankelAE_2017"="yellow3", "WindTT_2020" = "orchid1", "LeeKA_2022"="skyblue2","McCullochJA_2022"="black", "PetersBA_2019"="darkorange4", PernigoniN_2021="steelblue4"   ) ) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Country_sample), width=0.03,offset=0.1 ) + scale_fill_manual(values = c("GBR"="#6A3D9A" , "NLD"="#FB9A99" , "USA"="brown" , "CHE" = "#FDBF6F" ,  "Other"="grey" ) ) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1 ) + scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e")) +
  geom_balance(node=3760  , fill="grey", color=NA, alpha=0.3) -> Collinsella_melanoma #+
 # geom_balance(node=3757  , fill="skyblue", color=NA, alpha=0.3) #Clade + sister
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/Melanoma_collinsella.pdf",Collinsella_melanoma)
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/Melanoma_collinsella.png",Collinsella_melanoma)

#Make association with phenotypes
Info %>% filter(! is.na(nivolumab)) %>% filter(!is.na(Age)) %>% 
  glm(Melanoma_clade ~ nivolumab + antacids + H2RA + psychotropics + PPI + pembrolizumab + `peg-IFN` +methotrexate + `anti-PD1`, . , family=binomial()  ) %>% summary()
Info %>% filter(! is.na(nivolumab)) %>% filter(!is.na(Age)) %>% 
  glm(Melanoma_clade ~ previous_therapy , . , family=binomial()  ) %>% summary() #effect=1.9, P=0.089
Info %>% filter(! is.na(history_of_therapy)) %>% filter(!is.na(Age)) %>% 
  glm(Melanoma_clade ~ history_of_therapy , . , family=binomial()  ) %>% summary() #effect=-0.28, P=0.4

Info %>% glmer(melanoma ~ Melanoma_clade + Age + (1|Country) , . ,  family=binomial()) %>% summary()

Info %>% filter(melanoma == T & ! is.na(RECIST)) %>% 
  lme4::glmer( Melanoma_clade ~ RECIST + Age + (1|study_name), . , family=binomial() ) -> T1 
Info %>% filter(melanoma == T & ! is.na(RECIST)) %>% 
  lme4::glmer( Melanoma_clade ~  Age + (1|study_name), . , family=binomial() ) -> T0 
anova(T1, T0)
#Does association depend on abundance
#read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Analysis_melanoma/Tree_MAG/metaphlan/LeeKA_2022.tsv", skip=1) %>% filter(grepl("t__SGB14546_group",clade_name) ) -> Abundance_leeka
#Abundance_leeka %>% select(-clade_name) %>% gather(ID_anal, rel_abundance) %>% left_join(Info, . ) %>% filter(! is.na(rel_abundance))  %>% filter(melanoma == T & ! is.na(RECIST)) %>%
#  glm( as.factor(RECIST) ~ rel_abundance * Melanoma_clade + Age , . , family=binomial()  ) %>% summary()




Info %>% filter(Melanoma_clade == T) %>% group_by(melanoma) %>% summarise(n()) #123 Melanoma, 61 no melanoma
Info %>% filter(Melanoma_clade == T) %>% group_by(study_name, melanoma) %>% summarise(N = n()) %>% filter(melanoma==T)  #5 melanoma studies, 11 not melanoma
data.frame(c(123, 61), c(47, 3085)) %>% fisher.test()


Info %>% filter(Melanoma_clade == T) %>% filter(melanoma == F) %>% group_by(study_name) %>% summarise(n())
#PernigoniN_2021 is a study with prostate cancer
Info  %>% group_by(Melanoma_clade,study_name=="PernigoniN_2021") %>% summarise(N = n())  %>% spread(Melanoma_clade, N ) %>% as.data.frame() %>% column_to_rownames('study_name == "PernigoniN_2021"') %>% fisher.test()
Info %>% filter(,study_name=="PernigoniN_2021") %>% group_by(Melanoma_clade, Country) %>% summarise(n()) 
read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/PernigoniN_meta.tsv") -> Meta_per
Info %>% filter(study_name=="PernigoniN_2021") %>% select(ID_anal, Age, Country, Melanoma_clade) %>% left_join(Meta_per %>% rename(ID_anal = sample_id ) %>% select(-c(country, study_name)) ) -> Info_pernigoni
Info_pernigoni %>% glmer(Melanoma_clade ~ Age  + (1|location) + disease_subtype, . , family=binomial()  ) %>% summary()
Info_pernigoni %>% mutate(treatment = ifelse( is.na(treatment), "None", treatment ) )  %>% glm(Melanoma_clade ~ Age   + treatment, . , family=binomial()  ) %>% summary()
#any other cohort of interst?
Phenos %>% filter(ID_anal %in% filter(Info, Melanoma_clade==T)$ID_anal ) %>% group_by(disease, study_name) %>% summarise(n())
Info %>% left_join(., Phenos %>% select(ID_anal, disease) ) %>% filter(study_name=="LiJ_2014") %>% group_by(Melanoma_clade==T, disease ) %>% summarise(n()) 

#Melanoma, metabolic independence
Info$Melanoma_clade 
MetIndp = read_rds("/mnt/project/Make_Associations/Functional_enrichment/Metabolic_completness/MetabolicIndependence_data.rds")
MetIndp$ID_anal %>% sapply(function(x){ if ( grepl("__",x) ){ return(str_split(x, "__")[[1]][2])
  }else{ return(str_split(x, "_bin")[[1]][1]) }  } ) ->MetIndp$ID_anal
MetIndp$ID_anal %>% sapply(function(x){str_split(x, "_bin")[[1]][1] } ) ->MetIndp$ID_anal
MetIndp %>% left_join(Info) %>% filter(! is.na(Melanoma_clade)) -> MetIndp # %>% group_by(Melanoma_clade) %>% summarise(n())
Ass_indp = tibble()
for (M in unique(MetIndp$module)){
  MetIndp %>% filter(module == M) %>% lm(stepwise_module_is_complete ~  Melanoma_clade, .  ) %>% summary() -> DR
  DR$coefficients %>% as.data.frame() %>% rownames_to_column("Feature") %>% filter(Feature == "Melanoma_cladeTRUE") %>% mutate(Module = M) %>% rbind(Ass_indp, .) -> Ass_indp 
}
Ass_indp %>% as_tibble() %>% arrange(`Pr(>|t|)`) %>% mutate(FDR = p.adjust(`Pr(>|t|)`))
MetIndp %>% filter(module ==  "M00122") #Cobalamin : Check this study https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7073937/ : cobalamin levels in plasma and solid cancers
# A tibble: 2 × 2
#500FG has sample sin the melanoma clade
Info %>% filter(study_name == "SchirmerM_2016") %>% group_by(Melanoma_clade) %>% summarise(n())
#Melanoma_clade `n()`
#<lgl>          <int>
#1 FALSE            383
#2 TRUE              15
Analyse_500FG = function(File = "/mnt/project/Make_Associations/Phenotypes/Specific/Cytokines_500FG.tsv"){
  Info %>% select(-subject_id) %>% left_join( select(Phenos, c(ID_anal, subject_id) ) ) %>% filter(study_name == "SchirmerM_2016")  %>% select(Melanoma_clade, subject_id, Age, Sex) -> Analysis_500FG
  read_tsv(File) %>% rename(subject_id = ID_500FG) %>% left_join(Analysis_500FG, . ) -> Analysis_500FG
  Inverse_rank_normal = function(Measurement){qnorm((rank(Measurement,na.last="keep")-0.5)/sum(!is.na(Measurement))) }
  Results_cytokines = tibble()
  for (Cytokine in colnames(Analysis_500FG)[4:dim(Analysis_500FG)[2]] ){
    Analysis_500FG2 = Analysis_500FG 
    Analysis_500FG2[Cytokine] = Inverse_rank_normal(Analysis_500FG2[Cytokine])
    Formula = as.formula(paste0( "`", Cytokine, "` ~ Melanoma_clade + Age + Sex" ))
    Analysis_500FG2 %>% filter(! is.na( !!sym(Cytokine) ) )  %>% lm(Formula, . ) %>% summary() -> Res
    Res$coefficients["Melanoma_cladeTRUE",]  %>% t() %>% as_tibble() %>% mutate(Measurement = Cytokine) %>% rbind(Results_cytokines,.) -> Results_cytokines
  }
  Results_cytokines %>% mutate(FDR =  p.adjust(`Pr(>|t|)`, "fdr")  ) %>% arrange(`Pr(>|t|)`) -> res
  return(list(res, Analysis_500FG))
}
Analyse_500FG("/mnt/project/Make_Associations/Phenotypes/Specific/Cytokines_500FG.tsv")
Analyse_500FG("/mnt/project/Make_Associations/Phenotypes/Specific/CBcell.tsv")
Hormones_analysis = Analyse_500FG("/mnt/project/Make_Associations/Phenotypes/Specific/Hormones.tsv")
Hormones_analysis[[2]] %>% ggplot(aes(x=Melanoma_clade, y = B_TESC)) + geom_boxplot() + theme_bw() + ggforce::geom_sina()

#Select clades to compare
#1. All the rest
#2. Sister clade
Tree$tip.label[getDescendants(Tree, 3757)] -> clade_and_sister
clade_and_sister[!is.na(clade_and_sister)] -> clade_and_sister
Info %>% mutate(Melanoma_clade_sister = ifelse(ID_anal %in% Members_clade, F, ifelse(ID_anal %in% clade_and_sister, T ,F ) ) ) -> Info
#162 samples
Info %>% glm(Melanoma_clade_sister ~ Age + Country + melanoma, . , family=binomial() ) %>% summary()
#Check Gene presence/absence between clades
Prepare_genes = function(){
  #others
  DF1 = read_tsv("/mnt/project//Make_Associations/Functional_enrichment/MAG_analysis/Gene_presence/t__SGB14546_group.tsv")
  DF2 = read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Analysis_melanoma/Presence_gene_long.tsv")
  
  DF1 %>% select(Study, Sample, Gene ) %>% rename(ID_anal = Sample) -> DF1
  DF2 %>% rename(Study = Dataset) -> DF2
  rbind(DF1,DF2) -> DF
  DF %>% mutate(value = 1) %>% pivot_wider(names_from = Gene, values_from = value, values_fill = 0) %>% return()
  
}
Prepare_genes() -> Gene_info
Info %>% mutate(has_MAG = ifelse(ID_anal %in% Gene_info$ID_anal, T, F) ) -> Info
#Get prevalence in the clades of interest
Info %>% group_by(Melanoma_clade_sister,has_MAG) %>% summarise(n())
Info %>% group_by(Melanoma_clade, has_MAG) %>% summarise(n())

Info %>% filter(has_MAG==T, Melanoma_clade == T) %>% group_by(Country, study_name) %>% summarise(n())
#Country study_name       `n()`
#<chr>   <chr>            <int>
#  1 GBR     LeeKA_2022           5
#2 NLD     LeeKA_2022           3
#3 NLD     WindTT_2020          1
#4 USA     McCullochJA_2022     2
Info %>% filter(has_MAG==T) %>% group_by(Melanoma_clade) %>% summarise(n())
# A tibble: 2 × 2
#Melanoma_clade `n()`
#<lgl>          <int>
# 1 FALSE             50
#2 TRUE              11


Info %>% left_join(Gene_info) -> Info_with_genes
Info_with_genes %>% filter(has_MAG==T) -> Info_with_genes

Info_with_genes %>% filter(Melanoma_clade==T) %>% select(colnames(Gene_info)) %>% select(-c(Study, ID_anal)) %>% apply(2, function(x){ sum(x)/length(x) } ) -> Prevalences_melanoma
Info_with_genes %>% filter(Melanoma_clade==F) %>% select(colnames(Gene_info)) %>% select(-c(Study, ID_anal)) %>% apply(2, function(x){ sum(x)/length(x) } ) -> Prevalences_rest
abs(Prevalences_melanoma - Prevalences_rest) -> dif_Prevalence


Genes_test = names(dif_Prevalence)[dif_Prevalence>0.1]





N = 0
Results_test=tibble()
for (Gene in Genes_test ){
  #Add pseudocount
  Info_with_genes %>% filter(!is.na(!!sym(Gene))) -> Info_with_genes2
  mean(Info_with_genes2[Gene] %>% as_vector() ) -> Prevalence
  sum( filter(Info_with_genes2, Melanoma_clade==T  )[Gene] %>% as_vector()) -> Prevalence_high
  sum( filter(Info_with_genes2, Melanoma_clade==F  )[Gene] %>% as_vector()) -> Prevalence_low
  N = N + 1    
  cat(N ,Gene)
  Info_with_genes2 %>% group_by(!!sym(Gene), Melanoma_clade) %>% summarise(N = n()) %>% drop_na() %>% spread(Melanoma_clade, N) %>% ungroup() %>% as.data.frame() %>% column_to_rownames(Gene) -> Table_test
  if (dim(Table_test)[1] < 2){ next }
  Table_test[is.na(Table_test)] = 0
  Table_test = Table_test + 1 #Add pseudcount
  
  Test_p = fisher.test(Table_test)
  Prevalence = apply(Table_test, 2, function(x){ x[2]/sum(x)  }  )
  Fold_Change = as.vector(Prevalence[1])/as.vector(Prevalence[2])
  Results = tibble( Prevalence_rest = Prevalence[1], Prevalence_clade = Prevalence[2], Fold=Fold_Change, P = Test_p$p.value, Odds_ratio = Test_p$estimate, Lower_bound = Test_p$conf.int[1], Upper_bound = Test_p$conf.int[2]     ) 
  Results_test = rbind(Results_test  , Results %>% mutate(Gene=Gene) ) -> Results_test
}
Results_test %>% arrange(P) %>% mutate(FDR = p.adjust(P, "fdr") ) -> Results_test
write_tsv(Results_test, "/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Enrichment_melanomaClade.tsv")
Results_test %>% filter(FDR<0.05) -> Sign
Info_with_genes %>% select(ID_anal, Melanoma_clade, melanoma, Country) %>% as.data.frame() %>% column_to_rownames("ID_anal") %>% mutate(Melanoma_clade = as.numeric(Melanoma_clade), melanoma = as.numeric(melanoma) ) -> AnnotationColor
Info_with_genes  %>% select(ID_anal, Sign$Gene) %>% as.data.frame() %>% column_to_rownames("ID_anal") %>% t()  %>% pheatmap::pheatmap(show_colnames = F, fontsize_row = 4, annotation_col=AnnotationColor) 

#how many have annotation
#cazymes
Cazymes = read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Match_anpan/For_enrichment/Enrichment_tables/Cazymes.rds")
reaction_to_pathway <- as_tibble(stack(Cazymes))
as_tibble(reaction_to_pathway) %>% rename(Gene=values, Caz=ind) -> r_to_p
Results_test %>% left_join(r_to_p ) %>% group_by(is.na(Caz)) %>% summarise(n())
Results_test %>% distinct() %>% dplyr::mutate(Odds_ratio =  log10(Odds_ratio) ) %>%  dplyr::rename(`z value` = Odds_ratio )  -> Input5 
Run_enrichment(Input5,  Table=Cazymes )


#EC
read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/For_enrichment/UniRef90/ec.rds") -> ECs
ECs_tibble <- enframe(ECs, name = "EC", value = "UniRef90")
ECs_tibble %>% separate_rows(UniRef90, sep="\t") -> ECs_tibble
Enrichment_tables = read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/For_enrichment/metacyc_pathways.rds")
Translate_EC_to_reaction = read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis//For_enrichment/EC_to_reaction.tsv.gz")
Translate_EC_to_reaction %>%  separate_rows(Pathways, sep = ",") -> Translate_EC_to_reaction

Results_test %>% left_join(ECs_tibble %>% rename(Gene = UniRef90), relationship = "many-to-many" ) -> merged_ec
merged_ec %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #319/1,606 annotated

merged_ec %>% left_join(Translate_EC_to_reaction %>% rename(EC = ec)) %>% distinct() -> Merged_info

Run_enrichment(merged_ec %>% left_join(Translate_EC_to_reaction %>% rename(EC = ec)) %>% distinct() %>% dplyr::mutate(Odds_ratio =  log10(Odds_ratio) ) %>%  dplyr::rename(`z value` = Odds_ratio ) %>% mutate(Gene = Reaction)   ,  Table=Enrichment_tables ) -> R_reaction
R_reaction %>% as_tibble() %>% arrange(pval) %>% filter(padj < 0.05)


Enrichment_tables %>% stack() %>% as_tibble() %>% rename(Reaction =values , Pathways = ind  ) -> Pathwyas_tibble
Pathwyas_tibble %>% group_by(Pathways) %>% summarise(N = n()) -> Pathways_numbers
Results_test_pathway=tibble()
Completedness_total = tibble(ID_anal = Info_with_genes$ID_anal)
for (Path in unique(Merged_info$Pathways) ){
  print(Path)
  if(is.na(Path)){ next }
  Translate_EC_to_reaction %>% filter(Pathways == Path) -> ECs_to_check
  Denominator =  filter(Pathways_numbers, Pathways == Path)$N
  ECs_tibble %>% filter(EC %in% ECs_to_check$ec) -> ECs_to_check
  Info_with_genes  %>% select( c(ID_anal, one_of(ECs_to_check$UniRef90 ) ) ) ->N_pathway
  N_pathway %>% gather(UniRef90, value ,colnames(select(N_pathway, -ID_anal)) )  %>% left_join(ECs_to_check) %>% group_by(ID_anal, EC) %>% summarise(C = sum(value)) %>% mutate(Numerator = ifelse(C > 1 , 1, C ) ) %>% group_by(ID_anal) %>% summarise(N = sum(Numerator)) %>% 
    mutate(Fraction_complete = N/Denominator  ) %>% left_join( select(Info_with_genes, c(ID_anal,  Melanoma_clade) ) ) -> For_test
  wilcox.test( filter(For_test, Melanoma_clade == T)$Fraction_complete, filter(For_test, Melanoma_clade == F)$Fraction_complete) -> Result
  
  For_test %>% group_by(Melanoma_clade) %>% summarise(Median_completion = median(Fraction_complete) ) -> Result_info
  
  tibble(Pathway = Path, P=Result$p.value, Median_completion_clade = filter(Result_info,Melanoma_clade == T)$Median_completion,Median_completion_noclade = filter(Result_info,Melanoma_clade == F)$Median_completion   ) %>% rbind(Results_test_pathway, .) -> Results_test_pathway
  For_test %>% select(ID_anal, Fraction_complete) -> To_save
  colnames(To_save)[2] = Path
  left_join(Completedness_total, To_save, by="ID_anal") %>% as_tibble() -> Completedness_total
}
Results_test_pathway %>% arrange(P) %>% mutate( FDR = p.adjust(P, "fdr") ) -> Results_test_pathway
write_tsv(Results_test_pathway, "/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Enrichment_melanomaClade_completion.tsv")

Results_test_pathway %>% filter(Median_completion_clade>=0.75) %>% filter(FDR<0.05)
Results_test_pathway %>% filter(FDR<0.05) %>% filter(Median_completion_clade>Median_completion_noclade) %>% filter(Median_completion_clade >= 0.5 )


Results_test_pathway %>% filter(FDR<0.05) -> Sign
Translate_EC_to_reaction %>% filter(Pathways %in% Sign$Pathway) -> SignEC
ECs_tibble %>% filter(EC %in% SignEC$ec ) -> SignUniRef90

Info_with_genes  %>% select(ID_anal, Sign$Gene) %>% as.data.frame() %>% column_to_rownames("ID_anal") %>% t()  %>% pheatmap::pheatmap(show_colnames = F, fontsize_row = 4, annotation_col=AnnotationColor) 
Info_with_genes %>% select(ID_anal, one_of(SignUniRef90$UniRef90)) %>% as.data.frame() %>% column_to_rownames("ID_anal") %>% t()  %>% pheatmap::pheatmap(show_colnames = F, fontsize_row = 4, annotation_col=AnnotationColor) 
Completedness_total %>% as.data.frame() %>% column_to_rownames("ID_anal") %>% t()  %>% pheatmap::pheatmap(show_colnames = F, fontsize_row = 11, annotation_col=AnnotationColor) 

#KEGG
read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/For_enrichment/UniRef90/KO.rds") -> KO_table
KO_table %>% stack() %>% as_tibble() %>% rename(Gene = values, KO = ind) -> gene_to_ko
Results_test %>% left_join(gene_to_ko ) -> merged_ko
merged_ko %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #92/1,599 annotated

Run_enrichment(Input5  ,  Table=KO_table ) -> R_kegg
R_kegg  %>% arrange(pval) %>% filter(padj < 0.05)


#Go
read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/For_enrichment/UniRef90/Go.rds") -> go_terms
go_terms %>% stack() %>% as_tibble() %>% rename(Gene = values, Go = ind) -> gene_to_go
left_join(Results_test, gene_to_go) -> merged_go
merged_go %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #1.231

#enrichment of GO
Run_enrichment(Input5,  Table=go_terms ) -> R_go
R_go %>% as_tibble() %>% filter(padj < 0.05  )


library(clusterProfiler)
search_kegg_organism('caer', by='kegg_code')
kk <- enrichKEGG(keyType = "uniprot",
                  gene =  str_replace(filter(Results_test, Odds_ratio > 1, FDR<0.05)$Gene, "UniRef90_", "") ,
                 organism     = 'caer',
                 pvalueCutoff = 0.05, universe =str_replace(Results_test$Gene,"UniRef90_", "")   )

arrange(Results_test, desc(Odds_ratio)) -> N
N$Odds_ratio -> genelist_d
names(genelist_d) =  str_replace(N$Gene, "UniRef90_", "")

kk <- gseKEGG(keyType = "uniprot",
                 gene =  genelist_d,
                 organism     = 'caer',
                 pvalueCutoff = 0.05, verbose = T   )


head(kk)
xkk@result %>% as_tibble()
library(KEGGREST)
pathway_info <- keggGet("caer03010")





###########################################
##Functional analysis in centenarians#####
###########################################
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> Phenos
Phenos %>% select(ID_anal, study_name, IBD, Sex) %>% rename(sample_id = ID_anal) -> Phenos_merge


Process_RDS("/mnt/project/Make_Associations/Association/Results/Continent_stratified/Age,Country_Europe/Models/t__SGB4584/Age/Model.rds") -> Age_europe
Process_RDS("/mnt/project/Make_Associations/Association/Results/Continent_stratified/Age,Country_Asia/Models/t__SGB4584/Age/Model.rds") -> Age_Asia
Process_RDS("/mnt/project/Make_Associations/Association/Results/Cohort_specific/Age_XuQ_2021/Models/t__SGB4584/Age/Model.rds") -> Age_xu

rbind(Age_europe %>% mutate(Continent = "Europe") , Age_Asia %>% mutate(Continent = "Asia") ) %>% distinct(sample_id, .keep_all = T) -> Phylo_age
left_join(Phylo_age, Phenos_merge)  %>% distinct(sample_id, .keep_all = T) -> Phylo_age
#Read Tree
Tree_gnavus_all = Process_tree("/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB4584.TreeShrink.tre")
keep.tip(Tree_gnavus_all, Phylo_age$sample_id) -> Tree_gnavus_all
Phylo_age %>% mutate(Octagenarians = ifelse(Age>=80, T, F) ) -> Phylo_age
ggtree(Tree_gnavus_all, layout="fan", open.angle=15, size=0.1) %<+% Phylo_age  -> p
p + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Octagenarians), width=0.03,offset=0.1) + scale_fill_manual(values=c("TRUE"= "#E31A1C", "FALSE"="white") )  + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age), width=0.03,offset=0.1 )  + scale_fill_viridis_c(option = "viridis") + 
  #new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Country_sample), width=0.03,offset=0.1 ) + scale_fill_manual(values = c("GBR"="#6A3D9A" , "NLD"="#FB9A99" , "USA"="brown" , "Other"="grey" ) ) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1 ) + scale_fill_manual(values=c("Europe" = "#1f77b4", "Asia" = "#E31A1C")) +
  geom_balance(node=3982  , fill="grey", color=NA, alpha=0.3) -> Tree_elderClade

##simpler annotation tree
ggtree(Tree_gnavus_all, layout="fan", open.angle=15, size=0.1) %<+% Phylo_age  -> p
p + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age), width=0.03,offset=0.1 )  + scale_fill_viridis_c(option = "viridis") + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1 ) + scale_fill_manual(values=c("Europe" = "#1f77b4", "Asia" = "#E31A1C")) +
  geom_balance(node=3982  , fill="grey", color=NA, alpha=0.3) -> Tree_elderClade2


Tree_gnavus_all$tip.label[getDescendants(Tree_gnavus_all, 3982)] -> Aging_clade
Phylo_age %>% mutate(Elder_clade = ifelse(sample_id %in% Aging_clade, T, F ) ) -> Phylo_age
Phylo_age %>% mutate(Elder_clade = ifelse(Elder_clade == T, "Sample in clade", "Sample outside clade" ) ) %>% filter(Continent == "Europe") %>% filter(Age>18) %>% ggplot(aes(x=Elder_clade, y=Age)) + geom_boxplot() + ggforce::geom_sina(alpha=0.5) + theme_bw() + coord_flip() + xlab("Nonagenarian-enriched clade") + theme(text = element_text(size = 14), axis.text.y = element_text(angle = 90, hjust=0.5 ))  -> BoxPlotAgeEurope
Phylo_age %>% filter(Continent == "Europe") %>% glmer(Elder_clade ~ Age + (1|Country), . , family=binomial() ) %>% summary()
Phylo_age %>% filter(Continent == "Asia") %>% filter(Age < 90) %>% glmer(Elder_clade ~ Age + (1|Country), . , family=binomial() ) %>% summary()

Phylo_age %>% filter(study_name %in% c("LLD", "300OB") ) %>% select(sample_id, Elder_clade, Age, Sex) %>% print(n=100)

Phylo_age %>% filter(Elder_clade == T) %>% group_by(study_name) %>% summarise(n())

read_rds("/mnt/project/Make_Associations/Association/Results/Continent_stratified/Country_Asia/Models/t__SGB4584/Nonagenarian/Model.rds") -> gnavus_data
gnavus_data$model_input -> Info_gnavus
gnavus_data$pglmm_fit$summary() -> Fit_gnavus
print("Preparing information table")
Fit_gnavus %>% filter( grepl("phylo", variable))  %>% filter(! grepl("std", variable)) %>% mutate(sample_id = factor(c("-", as.character(Info_gnavus$sample_id)))) %>%
  select(sample_id, median) %>% rename(phylo_effect_median = median) -> Phylo
left_join(Info_gnavus, Phylo) %>% mutate(sample_id = as.character(sample_id)) -> Phylo
#add age and other phenos
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> Phenos
Phenos %>% select(study_name, ID_anal, Age) %>% rename(sample_id= ID_anal) %>% left_join(Phylo, . ) -> Info_gnavus
Info_gnavus %>% mutate(Age_group = ifelse(Age ==0, "<1", ifelse(Age < 15, "(0-15)", ifelse( Age < 30, "[15, 30)", ifelse( Age < 70, "[30, 70)", ">=70") )))) %>% mutate(Age_group = factor(Age_group, levels=c( "<1", "(0-15)", "[15, 30)", "[30, 70)", ">=70"  ))) -> Info_gnavus
Tree_gnavus = Process_tree("/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB4584.TreeShrink.tre")
keep.tip(Tree_gnavus, Info_gnavus$sample_id) -> Tree_gnavus
ggtree(Tree_gnavus, layout="fan", open.angle=15, size=0.1) %<+% Info_gnavus  -> p
getMRCA(Tree_gnavus, match( filter(Info_gnavus, Nonagenarian == T)$sample_id[20:40], Tree_gnavus$tip.label  ) )

p + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age), width=0.03,offset=0.1 )  + scale_fill_viridis_c(option = "viridis") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Nonagenarian), width=0.03,offset=0.1 )  + scale_fill_manual(values = c("FALSE"="#132b43", "TRUE"="#56b1f7")) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Country), width=0.03,offset=0.1 )  + #scale_fill_manual(values=c25)
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=study_name), width=0.03,offset=0.1 ) + scale_fill_manual(values=c25) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age_group), width=0.03,offset=0.1 ) + scale_fill_manual(values = c("<1" = "#5E1DB5",  "(0-15)"= "#5A3C82", "[15, 30)" = "#0E27E8", "[30, 70)"= "#EBB649",  ">=70"="#B5731D" ), na.value="white" ) +
  geom_balance(node=1874  , fill="grey", color=NA, alpha=0.3) -> Tree_elderClade_Asia
#simplify tree
ggtree(Tree_gnavus, layout="fan", open.angle=15, size=0.1) %<+% Info_gnavus  -> p
p + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age), width=0.03,offset=0.1 )  + scale_fill_viridis_c(option = "viridis") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Nonagenarian), width=0.03,offset=0.1 )  + scale_fill_manual(values = c("TRUE"="#d62728", "FALSE"="grey")) +
  geom_balance(node=1874  , fill="grey", color=NA, alpha=0.3) -> Tree_elderClade_Asia2



#Tree 2  
Info_gnavus2 = filter(Info_gnavus, study_name == "XuQ_2021" )
keep.tip(Tree_gnavus, Info_gnavus2$sample_id) -> Tree_gnavus2
Info_gnavus2 %>% mutate(Age_group = ifelse( Age < 60, "[50, 60)", ifelse( Age < 70, "[60, 70)" ,ifelse(Age<80, "[70, 80)", ifelse(Age<90, "[80, 90)", ">90" ) ) )))  -> Info_gnavus2

ggtree(Tree_gnavus2, layout="fan", open.angle=15, size=0.1) %<+% Info_gnavus2  -> p
p + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age), width=0.03,offset=0.1 )  + scale_fill_viridis_c(option = "viridis") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Nonagenarian), width=0.03,offset=0.1 )  + scale_fill_manual(values = c("FALSE"="#132b43", "TRUE"="#56b1f7")) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age_group), width=0.03,offset=0.1 ) + scale_fill_manual(values = c("[50, 60)" = "#5E1DB5",  "[60, 70)"= "#5A3C82", "[70, 80)" = "#0E27E8", "[80, 90)"= "#EBB649",  ">90"="#B5731D" ), na.value="white" ) + 
  geom_balance(node=137  , fill="grey", color=NA, alpha=0.3) 

Names = filter(Info_gnavus2, phylo_effect_median>1.5)$sample_id
match( Names , Tree_gnavus2$tip.label) -> CHECK
common_ancestor <- getMRCA(Tree_gnavus2, CHECK)

match( Names , Tree_gnavus$tip.label) -> CHECK
common_ancestor <- getMRCA(Tree_gnavus, CHECK)

match( Names , Tree_gnavus_all$tip.label) -> CHECK
common_ancestor <- getMRCA(Tree_gnavus_all, CHECK)


ggsave("/mnt/project/Make_Associations/Association/Results/Plots/Age_gnavus_ALL_s.pdf",Tree_elderClade2, width = 6, height = 6)
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/Age_gnavus_ALL_s.png",Tree_elderClade2)

ggsave("/mnt/project/Make_Associations/Association/Results/Plots/Age_gnavus_Asia_s.pdf",Tree_elderClade_Asia2, width = 12, height = 12)
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/Age_gnavus_Asia_s.png",Tree_elderClade_Asia2)

ggsave("/mnt/project/Make_Associations/Association/Results/Plots/Age_gnavus_Europe_boxplot.pdf",BoxPlotAgeEurope, width=3, height=3)
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/Age_gnavus_Europe_boxplot.png",BoxPlotAgeEurope, width=3, height=3)

Phylo_age %>% write_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/CladeRgnavus.tsv")
#Is there a depletion of IBD in this clade?
library(lmerTest)
Phylo_age %>% filter(Continent == "Europe")  -> European_IBD
European_IBD %>% mutate(IBD = ifelse(is.na(IBD), 0 , IBD ) ) -> European_IBD
European_IBD %>% glmer(Elder_clade ~ IBD + Age + (1|Country), . , family=binomial() ) %>% summary()
European_IBD %>% select(IBD, Elder_clade) %>% group_by(Elder_clade, IBD) %>% summarise(N = n()) %>% spread(IBD, N) %>% as.data.frame() %>% column_to_rownames("Elder_clade") %>% fisher.test()
#Clade     Yes No
#IBD Yes
#IBD No     

#Gene info R. gnavus
read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Gene_presence/t__SGB4584.tsv") -> UniRef90
UniRef90_wide <- UniRef90 %>% select(Sample, Gene) %>% filter(Sample %in% Phylo_age$sample_id) %>%
  mutate(N = 1) %>% distinct() %>% spread(Gene, N) 
UniRef90_wide[is.na(UniRef90_wide)] =  0

#Age_europe
#Age_Asia

#1. Check in European samples
Phylo_age %>% filter(sample_id %in% UniRef90_wide$Sample) %>% filter(Continent == "Europe") %>% left_join(., UniRef90_wide %>% rename(sample_id = Sample), by="sample_id" ) -> EuropeanUniRef
Phylo_age %>% distinct(sample_id, .keep_all=T)  -> Phylo_age

EuropeanUniRef %>% filter(Elder_clade==T) %>% select(colnames(UniRef90_wide%>% select(-Sample)))  %>% apply(2, function(x){ sum(x)/length(x) } ) -> Prevalences_cent
EuropeanUniRef %>% filter(Elder_clade==F) %>% select(colnames(UniRef90_wide %>% select(-Sample)))  %>% apply(2, function(x){ sum(x)/length(x) } ) -> Prevalences_nocent
abs(Prevalences_cent - Prevalences_nocent) -> dif_Prevalence_cent


Genes_test_c = names(dif_Prevalence_cent)[dif_Prevalence_cent>0.3]

N = 0
Results_test_c=tibble()
for (Gene in Genes_test_c ){
  #Add pseudocount
  EuropeanUniRef %>% filter(!is.na(!!sym(Gene))) -> Info_with_genes2
  mean(Info_with_genes2[Gene] %>% as_vector() ) -> Prevalence
  sum( filter(Info_with_genes2, Elder_clade==T  )[Gene] %>% as_vector()) -> Prevalence_high
  sum( filter(Info_with_genes2, Elder_clade==F  )[Gene] %>% as_vector()) -> Prevalence_low
  N = N + 1    
  cat(N ,Gene)
  Info_with_genes2 %>% group_by(!!sym(Gene), Elder_clade) %>% summarise(N = n()) %>% drop_na() %>% spread(Elder_clade, N) %>% ungroup() %>% as.data.frame() %>% column_to_rownames(Gene) -> Table_test
  if (dim(Table_test)[1] < 2){ next }
  Table_test[is.na(Table_test)] = 0
  Table_test = Table_test + 1 #Add pseudcount
  
  Test_p = fisher.test(Table_test)
  Prevalence = apply(Table_test, 2, function(x){ x[2]/sum(x)  }  )
  Fold_Change = as.vector(Prevalence[1])/as.vector(Prevalence[2])
  Results = tibble( Prevalence_rest = Prevalence[1], Prevalence_clade = Prevalence[2], Fold=Fold_Change, P = Test_p$p.value, Odds_ratio = Test_p$estimate, Lower_bound = Test_p$conf.int[1], Upper_bound = Test_p$conf.int[2]     ) 
  Results_test_c = rbind(Results_test_c  , Results %>% mutate(Gene=Gene) ) 
}
Results_test_c %>% arrange(P) %>% mutate(FDR = p.adjust(P, "fdr") ) -> Results_test_c
write_tsv(Results_test_c, "/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Enrichment_nonagenarianClade.tsv")
Results_test_c %>% filter(FDR<0.05) -> Sign
EuropeanUniRef %>% distinct(sample_id, .keep_all=T) %>% select(sample_id, Elder_clade, Age, Country) %>% as.data.frame() %>% column_to_rownames("sample_id") %>% mutate(Elder_clade = as.numeric(Elder_clade)) -> AnnotationColor
EuropeanUniRef %>% distinct(sample_id, .keep_all=T)  %>% select(sample_id, Sign$Gene) %>% as.data.frame() %>% column_to_rownames("sample_id") %>% t()  %>% pheatmap::pheatmap(show_colnames = F, fontsize_row = 4, annotation_col=AnnotationColor) 

#cazymes
Results_test_c %>% left_join(r_to_p ) %>% group_by(is.na(Caz)) %>% summarise(n())
Results_test_c %>% distinct() %>% dplyr::mutate(Odds_ratio =  log10(Odds_ratio) ) %>%  dplyr::rename(`z value` = Odds_ratio )  -> Input5 
Run_enrichment(Input5,  Table=Cazymes )
#EC
Results_test_c %>% left_join(ECs_tibble %>% rename(Gene = UniRef90), relationship = "many-to-many" ) -> merged_ec2
merged_ec2 %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #51/501 annotated
merged_ec2 %>% left_join(Translate_EC_to_reaction %>% rename(EC = ec)) %>% distinct() -> Merged_info2

Run_enrichment(merged_ec2 %>% left_join(Translate_EC_to_reaction %>% rename(EC = ec)) %>% distinct() %>% dplyr::mutate(Odds_ratio =  log10(Odds_ratio) ) %>%  dplyr::rename(`z value` = Odds_ratio ) %>% mutate(Gene = Reaction)   ,  Table=Enrichment_tables ) -> R_reaction2
R_reaction2 %>% as_tibble() %>% arrange(pval) %>% filter(padj < 0.05)

Results_test_pathway_c=tibble()
Completedness_total2 = tibble(sample_id = EuropeanUniRef$sample_id)
for (Path in unique(Merged_info2$Pathways) ){
  print(Path)
  if(is.na(Path)){ next }
  Translate_EC_to_reaction %>% filter(Pathways == Path) -> ECs_to_check
  Denominator =  filter(Pathways_numbers, Pathways == Path)$N
  ECs_tibble %>% filter(EC %in% ECs_to_check$ec) -> ECs_to_check
  EuropeanUniRef  %>% select( c(sample_id, one_of(ECs_to_check$UniRef90 ) ) ) ->N_pathway
  N_pathway %>% gather(UniRef90, value ,colnames(select(N_pathway, -sample_id)) )  %>% left_join(ECs_to_check) %>% group_by(sample_id, EC) %>% summarise(C = sum(value)) %>% mutate(Numerator = ifelse(C > 1 , 1, C ) ) %>% group_by(sample_id) %>% summarise(N = sum(Numerator)) %>% 
    mutate(Fraction_complete = N/Denominator  ) %>% left_join( select(EuropeanUniRef, c(sample_id,  Elder_clade) ) ) -> For_test
  wilcox.test( filter(For_test, Elder_clade == T)$Fraction_complete, filter(For_test, Elder_clade == F)$Fraction_complete) -> Result
  
  For_test %>% group_by(Elder_clade) %>% summarise(Median_completion = median(Fraction_complete) ) -> Result_info
  
  tibble(Pathway = Path, P=Result$p.value, Median_completion_clade = filter(Result_info,Elder_clade == T)$Median_completion,Median_completion_noclade = filter(Result_info,Elder_clade == F)$Median_completion   ) %>% rbind(Results_test_pathway_c, .) -> Results_test_pathway_c
  For_test %>% select(sample_id, Fraction_complete) -> To_save
  colnames(To_save)[2] = Path
  left_join(Completedness_total2, To_save, by="sample_id") %>% as_tibble() -> Completedness_total2
}
Results_test_pathway_c %>% arrange(P) %>% mutate( FDR = p.adjust(P, "fdr") ) -> Results_test_pathway_c
write_tsv(Results_test_pathway, "/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Enrichment_elderClade_completionEurope.tsv")

Results_test_pathway %>% filter(Median_completion_clade>=0.75) %>% filter(FDR<0.05)
Results_test_pathway %>% filter(FDR<0.05) %>% filter(Median_completion_clade>Median_completion_noclade) %>% filter(Median_completion_clade >= 0.5 )


Results_test_pathway %>% filter(FDR<0.05) -> Sign
Translate_EC_to_reaction %>% filter(Pathways %in% Sign$Pathway) -> SignEC
ECs_tibble %>% filter(EC %in% SignEC$ec ) -> SignUniRef90

Info_with_genes  %>% select(ID_anal, Sign$Gene) %>% as.data.frame() %>% column_to_rownames("ID_anal") %>% t()  %>% pheatmap::pheatmap(show_colnames = F, fontsize_row = 4, annotation_col=AnnotationColor) 
Info_with_genes %>% select(ID_anal, one_of(SignUniRef90$UniRef90)) %>% as.data.frame() %>% column_to_rownames("ID_anal") %>% t()  %>% pheatmap::pheatmap(show_colnames = F, fontsize_row = 4, annotation_col=AnnotationColor) 
Completedness_total %>% as.data.frame() %>% column_to_rownames("ID_anal") %>% t()  %>% pheatmap::pheatmap(show_colnames = F, fontsize_row = 11, annotation_col=AnnotationColor) 

#KEGG
Results_test_c %>% left_join(gene_to_ko ) -> merged_ko
merged_ko %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #4 annotated
Run_enrichment(Input5  ,  Table=KO_table ) -> R_kegg
R_kegg  %>% arrange(pval) %>% filter(padj < 0.05)
#Go
left_join(Results_test_c, gene_to_go) -> merged_go
merged_go %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #1.231
#enrichment of GO
Run_enrichment(Input5,  Table=go_terms ) -> R_go
R_go %>% as_tibble() %>% filter(padj < 0.05  )



###Analysis XuQ
#Blast iso-
BlastResult = read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Analysis_nonagenarianClade/Blast_result.tsv")
Phylo_age %>% filter(sample_id %in% BlastResult$Sample) %>% left_join(BlastResult %>% rename(sample_id =Sample)) %>% drop_na() %>% select(-c(Identity, Length,Evalue)) %>%
  mutate(Presence=1) %>% spread(Protein, Presence) -> PhyloAgeBlast
#Present in all samples
#UniRef90

#gutSmash
GS_qin = read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/GutSmash/GutSmash_QinN_2021.tsv", col_names=F)
GS_qin %>% rename(sample_id = X1) %>% select(-X2) %>% distinct() %>% mutate(Presence=1) %>% spread(X3, Presence) -> GS_qin
GS_qin[is.na(GS_qin)] = 0
GS_qin %>% select(-sample_id) %>% apply(2, function(x){ sum(x)/length(x) } )
left_join(Phylo_age, GS_qin) %>% drop_na() -> Info_with_gutsmash
Results_gutsmash = tibble()
for (Module in colnames(GS_qin %>% select(-sample_id))){
  #Add pseudocount
  Info_with_gutsmash %>% filter(!is.na(!!sym(Module))) -> Info_with_genes2
  mean(Info_with_genes2[Module] %>% as_vector() ) -> Prevalence
  sum( filter(Info_with_genes2, Elder_clade==T  )[Module] %>% as_vector()) -> Prevalence_high
  sum( filter(Info_with_genes2, Elder_clade==F  )[Module] %>% as_vector()) -> Prevalence_low
  Info_with_genes2 %>% group_by(!!sym(Module), Elder_clade) %>% summarise(N = n()) %>% drop_na() %>% spread(Elder_clade, N) %>% ungroup() %>% as.data.frame() %>% column_to_rownames(Module) -> Table_test
  if (dim(Table_test)[1] < 2){ next }
  Table_test[is.na(Table_test)] = 0
  Table_test = Table_test + 1 #Add pseudcount
  
  Test_p = fisher.test(Table_test)
  Prevalence = apply(Table_test, 2, function(x){ x[2]/sum(x)  }  )
  Fold_Change = as.vector(Prevalence[1])/as.vector(Prevalence[2])
  Results = tibble( Prevalence_rest = Prevalence[1], Prevalence_clade = Prevalence[2], Fold=Fold_Change, P = Test_p$p.value, Odds_ratio = Test_p$estimate, Lower_bound = Test_p$conf.int[1], Upper_bound = Test_p$conf.int[2]     ) 
  Results_gutsmash = rbind(Results_gutsmash  , Results %>% mutate(BGC=Module) ) 
}
Results_gutsmash %>% as_tibble() %>% arrange(P) %>% mutate(FDR= p.adjust(P, "fdr") )

#Replicate in rest of non-Qin samples
GS_no_qin = read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/GutSmash/GutSmash_noQinN_2021.tsv", col_names=F)
GS_no_qin %>% rename(sample_id = X1) %>% select(-X2) %>% distinct() %>% mutate(Presence=1) %>% spread(X3, Presence) -> GS_no_qin
GS_no_qin[is.na(GS_no_qin)] = 0
left_join(Phylo_age, GS_no_qin  ) %>% drop_na() -> Info_with_gutsmash2
Info_with_gutsmash2 %>% group_by(Elder_clade) %>% summarise(n())
Info_with_gutsmash2 %>% group_by(!!sym("BGC type: Flavoenzyme_sugar_catabolism"), Elder_clade) %>% summarise(N = n()) %>% drop_na() %>% spread(Elder_clade, N) %>% ungroup() %>% as.data.frame() %>% column_to_rownames(Module) -> Table_test
Table_test[is.na(Table_test)] = 0 #; Table_test = Table_test + 1
Test_p = fisher.test(Table_test)
print(Test_p$p.value) ; print(Test_p$estimate)


#UniRef90
read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Analysis_nonagenarianClade/Presence_gene_long.tsv") -> UniRef90_xu
UniRef90_wide_xu <- UniRef90_xu %>% select(ID_anal, Gene) %>% filter(ID_anal %in% Phylo_age$sample_id) %>%
  mutate(N = 1) %>% distinct() %>% spread(Gene, N) 
UniRef90_wide_xu[is.na(UniRef90_wide_xu)] =  0
Phylo_age %>% filter(sample_id %in% UniRef90_wide_xu$ID_anal) -> XuUniRef
Phylo_age %>% distinct(sample_id, .keep_all=T)  -> Phylo_age

XuUniRef %>% left_join(UniRef90_wide_xu %>% rename(sample_id = ID_anal)) -> XuUniRef

XuUniRef %>% filter(Elder_clade==T) %>% select(colnames(UniRef90_wide_xu%>% select(-ID_anal)))  %>% apply(2, function(x){ sum(x)/length(x) } ) -> Prevalences_cent
XuUniRef %>% filter(Elder_clade==F) %>% select(colnames(UniRef90_wide_xu %>% select(-ID_anal)))  %>% apply(2, function(x){ sum(x)/length(x) } ) -> Prevalences_nocent
abs(Prevalences_cent - Prevalences_nocent) -> dif_Prevalence_cent
Genes_test_c = names(dif_Prevalence_cent)[dif_Prevalence_cent>0.3]
N = 0
Results_test_xu=tibble()
for (Gene in Genes_test_c ){
  #Add pseudocount
  XuUniRef %>% filter(!is.na(!!sym(Gene))) -> Info_with_genes2
  mean(Info_with_genes2[Gene] %>% as_vector() ) -> Prevalence
  sum( filter(Info_with_genes2, Elder_clade==T  )[Gene] %>% as_vector()) -> Prevalence_high
  sum( filter(Info_with_genes2, Elder_clade==F  )[Gene] %>% as_vector()) -> Prevalence_low
  N = N + 1    
  cat(N ,Gene)
  Info_with_genes2 %>% group_by(!!sym(Gene), Elder_clade) %>% summarise(N = n()) %>% drop_na() %>% spread(Elder_clade, N) %>% ungroup() %>% as.data.frame() %>% column_to_rownames(Gene) -> Table_test
  if (dim(Table_test)[1] < 2){ next }
  Table_test[is.na(Table_test)] = 0
  Table_test = Table_test + 1 #Add pseudcount
  
  Test_p = fisher.test(Table_test)
  Prevalence = apply(Table_test, 2, function(x){ x[2]/sum(x)  }  )
  Fold_Change = as.vector(Prevalence[1])/as.vector(Prevalence[2])
  Results = tibble( Prevalence_rest = Prevalence[1], Prevalence_clade = Prevalence[2], Fold=Fold_Change, P = Test_p$p.value, Odds_ratio = Test_p$estimate, Lower_bound = Test_p$conf.int[1], Upper_bound = Test_p$conf.int[2]     ) 
  Results_test_xu = rbind(Results_test_xu  , Results %>% mutate(Gene=Gene) ) 
}
Results_test_xu %>% arrange(P) %>% mutate(FDR = p.adjust(P, "fdr") ) -> Results_test_xu
write_tsv(Results_test_xu, "/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Enrichment_nonagenarianClade_xu.tsv")

Results_test_xu %>% mutate(Cohorts = "XuQ_2021") %>% rbind(Results_test_c %>% mutate(Cohorts = "European")) -> AllEnrichment
write_tsv(AllEnrichment, "/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Enrichment_elderClade_XuANDEurope.tsv")



Results_test_xu %>% filter(FDR<0.05) -> Sign
XuUniRef %>% distinct(sample_id, .keep_all=T) %>% select(sample_id, Elder_clade, Age, Country) %>% as.data.frame() %>% column_to_rownames("sample_id") %>% mutate(Elder_clade = as.numeric(Elder_clade)) -> AnnotationColor
XuUniRef %>% distinct(sample_id, .keep_all=T)  %>% select(sample_id, Sign$Gene) %>% as.data.frame() %>% column_to_rownames("sample_id") %>% t()  %>% pheatmap::pheatmap(show_colnames = F, fontsize_row = 4, annotation_col=AnnotationColor) 

read_tsv( "/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Enrichment_nonagenarianClade.tsv") -> European
European %>% filter(FDR<0.05) %>% group_by( Gene %in% Sign$Gene)  %>% summarise(n())
Sign %>% filter(FDR<0.05) %>%  group_by( Gene %in% filter(European, FDR<0.05)$Gene ) %>% summarise(n()) #287 replicated in European clade


#EC
Results_test_xu %>% left_join(ECs_tibble %>% rename(Gene = UniRef90), relationship = "many-to-many" ) -> merged_ecxu
merged_ecxu %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #82/821 annotated
merged_ecxu %>% left_join(Translate_EC_to_reaction %>% rename(EC = ec)) %>% distinct() -> Merged_info3

Run_enrichment(merged_ecxu %>% left_join(Translate_EC_to_reaction %>% rename(EC = ec)) %>% distinct() %>% dplyr::mutate(Odds_ratio =  log10(Odds_ratio) ) %>%  dplyr::rename(`z value` = Odds_ratio ) %>% mutate(Gene = Reaction)   ,  Table=Enrichment_tables ) -> R_reaction2
R_reaction2 %>% as_tibble() %>% arrange(pval) %>% filter(padj < 0.05)

Results_test_pathway_xu=tibble()
Completedness_total_xu = tibble(sample_id = XuUniRef$sample_id)
for (Path in unique(Merged_info3$Pathways) ){
  print(Path)
  if(is.na(Path)){ next }
  Translate_EC_to_reaction %>% filter(Pathways == Path) -> ECs_to_check
  Denominator =  filter(Pathways_numbers, Pathways == Path)$N
  ECs_tibble %>% filter(EC %in% ECs_to_check$ec) -> ECs_to_check
  XuUniRef  %>% select( c(sample_id, one_of(ECs_to_check$UniRef90 ) ) ) ->N_pathway
  N_pathway %>% gather(UniRef90, value ,colnames(select(N_pathway, -sample_id)) )  %>% left_join(ECs_to_check) %>% group_by(sample_id, EC) %>% summarise(C = sum(value)) %>% mutate(Numerator = ifelse(C > 1 , 1, C ) ) %>% group_by(sample_id) %>% summarise(N = sum(Numerator)) %>% 
    mutate(Fraction_complete = N/Denominator  ) %>% left_join( select(XuUniRef, c(sample_id,  Elder_clade) ) ) -> For_test
  wilcox.test( filter(For_test, Elder_clade == T)$Fraction_complete, filter(For_test, Elder_clade == F)$Fraction_complete) -> Result
  
  For_test %>% group_by(Elder_clade) %>% summarise(Median_completion = median(Fraction_complete) ) -> Result_info
  
  tibble(Pathway = Path, P=Result$p.value, Median_completion_clade = filter(Result_info,Elder_clade == T)$Median_completion,Median_completion_noclade = filter(Result_info,Elder_clade == F)$Median_completion   ) %>% rbind(Results_test_pathway_xu, .) -> Results_test_pathway_xu
  For_test %>% select(sample_id, Fraction_complete) -> To_save
  colnames(To_save)[2] = Path
  left_join(Completedness_total_xu, To_save, by="sample_id") %>% as_tibble() -> Completedness_total_xu
}
Results_test_pathway_xu %>% arrange(P) %>% mutate( FDR = p.adjust(P, "fdr") ) -> Results_test_pathway_xu
write_tsv(Results_test_pathway_xu, "/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Enrichment_elderClade_completionXu.tsv")
Results_test_pathway_xu %>% filter(FDR<0.05) %>% filter(Pathway %in% filter(Results_test_pathway_c, FDR<0.05 )$Pathway  )

Results_test_pathway_xu %>% mutate(Cohorts = "XuQ_2021") %>% rbind(Results_test_pathway %>% mutate(Cohorts = "European")) -> AllCompletedness
write_tsv(AllCompletedness, "/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Enrichment_elderClade_completion_XuANDEurope.tsv")


Completedness_total_xu %>% left_join(Phylo_age) %>% mutate(Nonagenarian = ifelse(Age<=90, T, F) ) %>% ggplot(aes(y=`PWY0-1315`, x=Elder_clade, col=Nonagenarian )) + geom_boxplot() + theme_bw() + geom_jitter()
Completedness_total2 %>% left_join(Phylo_age) %>% mutate(Nonagenarian = ifelse(Age<=90, T, F) ) %>% ggplot(aes(y=`PWY0-1315`, x=Elder_clade, col=Nonagenarian )) + geom_boxplot() + theme_bw() + geom_jitter()


#Go
left_join(Results_test_xu, gene_to_go) -> merged_go
merged_go %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #562/821
#enrichment of GO
Run_enrichment(Input5,  Table=go_terms ) -> R_go
R_go %>% as_tibble() %>% filter(padj < 0.05  )





#########################################
#########Functional analysis on BMI#####
########################################
read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Gene_presence/t__SGB4285_group.tsv") -> UniRef90
read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Gene_presence/ec/t__SGB4285_group.tsv") -> EC

read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Bromii_BMI/BMI_Ruminococcus_bromii.tsv") -> Bromii_BMI

UniRef90 %>% filter(Sample %in% Bromii_BMI$sample_id ) -> UniRef90
Bromii_BMI %>% mutate(MAG_available = ifelse(sample_id %in% UniRef90$Sample, T, F) ) -> Bromii_BMI

Tree_bromii = Process_tree("/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB4285_group.TreeShrink.tre")
keep.tip(Tree_bromii,Bromii_BMI$sample_id) -> Tree_bromii

Bromii_BMI %>% mutate(BMI_group = ifelse(BMI<18.5, "BMI<18.5", ifelse(BMI>=18.5 & BMI<25, "18.5≤BMI<25", ifelse(BMI>=25 & BMI<30, "25≤BMI<30", "30≤BMI" ) ) ) ) -> Bromii_BMI
Bromii_BMI = Bromii_BMI %>% mutate(Age_group = ifelse(Age ==0, "<1", ifelse(Age < 15, "(0-15)", ifelse( Age < 30, "[15, 30)", ifelse( Age < 70, "[30, 70)", ">=70") )))) %>% mutate(Age_group = factor(Age_group, levels=c( "<1", "(0-15)", "[15, 30)", "[30, 70)", ">=70"  )))



#show tree only on a subset
set.seed(12)
sample(Bromii_BMI$sample_id, 1000) -> SubSet
Bromii_BMI %>% filter(sample_id %in% SubSet) -> Bromii_BMI2
keep.tip(Tree_bromii, SubSet) -> Tree_bromii2

ggtree(Tree_bromii2, layout = "fan",  open.angle=15, size=0.1 ) %<+% Bromii_BMI2 -> p #+
  p + geom_fruit( geom="geom_tile", mapping = aes(fill=BMI_group), width=0.03,offset=0.1  ) + scale_fill_manual(values = c("BMI<18.5" = "blue",  "18.5≤BMI<25"= "grey", "25≤BMI<30" = "#ff7f0e", "30≤BMI"= "red"), na.value="white" ) +   labs(fill = "BMI") + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age_group), width=0.03,offset=0.1  ) + scale_fill_manual(values = c("<1" = "#5E1DB5",  "(0-15)"= "#5A3C82", "[15, 30)" = "#0E27E8", "[30, 70)"= "#EBB649",  ">=70"="#B5731D" ), na.value="white" ) +   labs(fill = "Age") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
 geom_balance(node=1003  , fill="skyblue", color=NA, alpha=0.3) + geom_balance(node=1683  , fill="#FF7F7F", color=NA, alpha=0.3) + geom_balance(node=1577  , fill="grey", color=NA, alpha=0.3)+ geom_balance(node=1577  , fill="grey", color=NA, alpha=0.3) + 
    #geom_balance(node=1573, fill="pink", color=NA, alpha=0.3) +
    new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=MAG_available), width=0.03,offset=0.1  ) + scale_fill_manual(values = c("FALSE" = "grey",  "TRUE"= "black"), na.value="white"  ) + labs(fill = "MAG") 
    
    
    
getMRCA(Tree_bromii2, seq(790, 805))

Find_node_ancestor = function(Phylum_name = "Firmicutes", subsample=50, SGB_tree=SGB_tree2){
  if (length(leafs_ids) > subsample){ sample(leafs_ids, subsample) -> leafs_ids } 
  leaf_nodes <- match(leafs_ids, Tree_bromii2$tip.label)
  #thresh = mean(leaf_nodes) + 3*sd(leaf_nodes)
  #leaf_nodes = leaf_nodes[ leaf_nodes < thresh   ]

  common_ancestor <- getMRCA(Tree_bromii, leaf_nodes)
  print(leaf_nodes) ; print(common_ancestor)
  tibble(node =common_ancestor  ) %>% return()
}

Bromii_BMI2 %>% arrange(phylo_effect_median) %>% head(n=10)  %>% select(sample_id) %>% as_vector() ->   leafs_ids   
leaf_nodes <- match(leafs_ids, Tree_bromii$tip.label)
common_ancestor1 <- getMRCA(Tree_bromii, leaf_nodes)

Bromii_BMI2 %>% arrange(desc(phylo_effect_median)) %>% head(n=10)  %>% select(sample_id) %>% as_vector() ->   leafs_ids   
leaf_nodes <- match(leafs_ids, Tree_bromii$tip.label)
common_ancestor2 <- getMRCA(Tree_bromii, leaf_nodes)

common_ancestor3 <- getMRCA(Tree_bromii, seq(1578, 1590))

common_ancestor4 <- getMRCA(Tree_bromii, seq(1578, 1590))
common_ancestor5 <- getMRCA(Tree_bromii, c(449, 542))



ggtree(Tree_bromii, layout = "fan",  open.angle=15, size=0.1 ) %<+% Bromii_BMI -> p #+
p + geom_fruit( geom="geom_tile", mapping = aes(fill=BMI_group), width=0.03,offset=0.1  ) + scale_fill_manual(values = c("BMI<18.5" = "blue",  "18.5≤BMI<25"= "grey", "25≤BMI<30" = "#ff7f0e", "30≤BMI"= "red"), na.value="white" ) +   labs(fill = "BMI") + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age_group), width=0.03,offset=0.1  ) + scale_fill_manual(values = c("<1" = "#5E1DB5",  "(0-15)"= "#5A3C82", "[15, 30)" = "#0E27E8", "[30, 70)"= "#EBB649",  ">=70"="#B5731D" ), na.value="white" ) +   labs(fill = "Age") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  geom_balance(node=4476  , fill="skyblue", color=NA, alpha=0.3) + geom_balance(node=7517  , fill="#FF7F7F", color=NA, alpha=0.3) + geom_balance(node=7517  , fill="grey", color=NA, alpha=0.3)  + 
  geom_balance(node=6138, fill="pink", color=NA, alpha=0.3) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=MAG_available), width=0.03,offset=0.1  ) + scale_fill_manual(values = c("FALSE" = "grey",  "TRUE"= "black"), na.value="white"  ) + labs(fill = "MAG") -> ruminoBromi_tree


ggsave("/mnt/project/Make_Associations/Association/Results/Plots/BMI_bromi.pdf",ruminoBromi_tree)
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/BMI_bromi.png",ruminoBromi_tree)


Tree_bromii$tip.label[getDescendants(Tree_bromii, node = common_ancestor1)] -> Low_BMI
Tree_bromii$tip.label[getDescendants(Tree_bromii, node = common_ancestor2)] -> High_BMI
Tree_bromii$tip.label[getDescendants(Tree_bromii, node = common_ancestor5)] -> Low2_BMI



Bromii_BMI %>% mutate( Clade = ifelse(sample_id %in% Low_BMI, "Low", ifelse(sample_id %in% High_BMI, "High", ifelse(sample_id %in% Low2_BMI, "Low2", "Other" ) )  )   ) %>% drop_na() -> ToTest
Phenotypes %>% select(ID_anal, study_name, Sex) %>% rename(sample_id = ID_anal ) %>% left_join(ToTest, . ) -> ToTest

ToTest %>% mutate(study_name = ifelse(study_name == "LLD2", "LLD", study_name) ) -> ToTest
ToTest %>% mutate( study_name = ifelse(study_name=="300OB", "KurilshikovA_2019", ifelse(study_name=="DAG3", "GacesaR_2022", ifelse(study_name=="LLD", "ZhernakovaA_2016", ifelse(study_name=="IBD", "VilaAV_2018", study_name ) ) ) ) ) -> ToTest


ToTest %>% filter(Country %in% c("DNK", "GBR", "NLD") ) %>%  ggplot(aes(x=Clade, y=BMI)) + geom_boxplot(outlier.shape = NA, aes(fill=Clade)) + theme_bw() + geom_sina(alpha=0.5 ) + facet_wrap(~Country , scales="free") +  
  scale_fill_manual(values = c("High"= "#FF7F7F", "Low" = "skyblue", "Low2"="pink", "Other"="grey" ) ) + theme(text = element_text(size = 14), axis.text.y = element_blank()) + coord_flip() + theme(text= element_text(size = 14)) -> boxplots_bmi
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/BMI_bromi_boxplot.pdf",boxplots_bmi)
ggsave("/mnt/project/Make_Associations/Association/Results/Plots/BMI_bromi_boxplot.png",boxplots_bmi)

ToTest %>% filter(Clade != "Other")  %>% rbind(. , ToTest %>% mutate(study_name = "Overall") ) %>% mutate(study_name = factor(study_name, c("Overall", unique(ToTest$study_name)  ) ) ) %>% filter(Country %in% c("DNK", "GBR", "NLD") ) %>%  ggplot(aes(x=Clade, y=BMI, col=study_name  )) + geom_half_boxplot(outlier.shape = NA) + theme_bw() + geom_sina(alpha=0.5 ) + facet_wrap(~Country , scales="free") +  
  theme(text = element_text(size = 14) ) + coord_flip() + scale_color_manual(values = c(c25[1:13], "brown" ) )

ToTest  %>% rbind(. , ToTest %>% mutate(study_name = "zOverall") ) %>% filter(Country %in% c("DNK", "GBR", "NLD") ) %>%  ggplot(aes(x=Clade, y=BMI, fill=study_name  )) + geom_half_boxplot(outlier.shape = NA) + theme_bw() +  facet_wrap(~Country , scales="free") +  geom_dotplot(binaxis = "y", method="histodot", stackdir="up", position = PositionDodge, dotsize = 0.1) + 
  +   theme(text = element_text(size = 14) )  + scale_fill_manual(values =c25) + coord_flip()



ToTest %>% filter(Country %in% c("DNK", "GBR", "NLD") ) %>% filter(! study_name %in% c("KarlsoonFH_2013", "LeeKA_2022") )  %>% drop_na() %>% lmerTest::lmer( BMI ~ Clade + Age + Sex + (1|Country) + (1|Country:study_name), . ) %>% summary()
ToTest  %>% drop_na() %>% lmerTest::lmer( BMI ~ Clade + Age + Sex + (1|Country) + (1|Country:study_name), . ) %>% summary()


ToTest %>% group_by(Clade, MAG_available) %>% summarise(n())

###I can do the enrichment at the Gene level or at the product level. Gene level is way bigger.
Do_enrichment_uniref90 = function( Feature = "Gene", Subset = T, Uniref = UniRef90, Test ="logistic" ){
  #Feature = "Product"
  if (Subset == T) { Uniref = Uniref %>% filter(Sample %in%   filter(ToTest, Clade!="Other")$sample_id ) ->  Uniref }
  Uniref %>% filter(Sample %in% ToTest$sample_id  )  %>% select( c("Sample", Feature)) %>% distinct(Sample, !!sym(Feature) , .keep_all = T) %>% mutate(Presence = 1) %>% spread( Feature, "Presence") -> UnirefTest
  #All missing values should be 0 for testing
  UnirefTest[is.na(UnirefTest)] = 0 
  #Add sample information
  UnirefTest %>% left_join( . ,Uniref %>% distinct(Sample, .keep_all =T ) %>% select(Study, Sample, Completeness, Contamination) )  -> UnirefTest
  UnirefTest %>% mutate(Clade = ifelse(Sample %in% Low2_BMI, "Low2", ifelse(Sample %in% High_BMI, "High", "Other"  ) ) ) -> UnirefTest
  #Keep only one entry per sample
  UnirefTest %>% distinct(Sample, .keep_all =T) -> UnirefTest
  UnirefTest %>% select(-c("Study", "Sample", "Completeness", "Contamination", "Clade")) %>% apply(2,function(x){ mean(x) } ) -> Prevalences
  

  Results_test = tibble()
  N = 0
  for (Gene in names(Prevalences[Prevalences < 0.9 & Prevalences > 0.1]) ){
    #Add pseudocount
    #if (Gene %in% c("Study", "Sample", "Completeness", "Contamination", "Clade") ){ next }
    mean(UnirefTest[Gene] %>% as_vector() ) -> Prevalence
    sum( filter(UnirefTest, Clade=="Low2"  )[Gene] %>% as_vector()) -> Prevalence_low
    sum( filter(UnirefTest, Clade=="High"  )[Gene] %>% as_vector()) -> Prevalence_high
    sum( filter(UnirefTest, Clade=="Other"  )[Gene] %>% as_vector()) -> Prevalence_other
    
    #if (Prevalence <= 0.1 | Prevalence >= 0.9 ){ next }
    N = N + 1    
    print( c(N ,Gene) )
    
    if (Test == "Logistic"){
    
    if (Prevalence_high == 0  ){  tibble(Feature = c("CladeLow2", "CladeOther") , Estimate=NA, `Std. Error`=NA, `z value`=NA,  `Pr(>|z|)`=0, Gene=Gene, Prevalence_low = Prevalence_low, Prevalence_high = Prevalence_high, Prevalence_other=Prevalence_other ) %>% rbind(Results_test) ->Results_test ;  next() }
    Formula = paste0("`",Gene, "` ~ Clade + Completeness + Contamination") 
    UnirefTest %>% glm(Formula , . , family=binomial()) %>% summary() -> Results
    Results$coefficients %>% as.data.frame() %>% rownames_to_column("Feature")  %>% filter(Feature %in%  c("CladeLow2", "CladeOther" )) %>% as_tibble() %>% mutate(Gene = Gene, Prevalence_low = Prevalence_low, Prevalence_high = Prevalence_high, Prevalence_other=Prevalence_other ) %>% rbind(Results_test) -> Results_test
    } else {
      UnirefTest %>% group_by(!!sym(Gene), Clade) %>% summarise(N = n())  %>% spread(Clade, N) %>% as.data.frame() %>% column_to_rownames(Gene) -> Table_test
      Table_test[is.na(Table_test)] = 0
    
      Table_test = Table_test + 1 #Add pseudcount
      for (Comparison in c("Other", "Low2") ){
        Table_test %>% select("High", Comparison) -> For_Test
        Test_p = fisher.test(For_Test)
        Prevalence = apply(For_Test, 2, function(x){ x[2]/sum(x)  }  )
        Fold_Change = as.vector(Prevalence[1])/as.vector(Prevalence[2])
        Results = tibble( Comparison = Comparison, Prevalence_high = Prevalence[1], Prevalence_comparison = Prevalence[2], Fold=Fold_Change, P = Test_p$p.value, Odds_ratio = Test_p$estimate, Lower_bound = Test_p$conf.int[1], Upper_bound = Test_p$conf.int[2]     ) 
        Results_test = rbind(Results_test  , Results %>% mutate(Gene=Gene) ) -> Results_test
      }
    }
  }
  if ( Test == "Logistic" ){ Results_test  %>% mutate(FDR = p.adjust(`Pr(>|z|)`, "fdr" ) ) -> Results_test
  } else { Results_test  %>% mutate(FDR = p.adjust(P, "fdr" ) ) -> Results_test }
  
  #Make Heatmap
  UnirefTest %>% select(Sample, Clade ) %>% as.data.frame()   %>% column_to_rownames("Sample") -> AnnotationColor
  Sign = filter(Results_test, FDR<0.05)$Gene
  
  if (Subset != T){
    #1. With 'Other' clade
    UnirefTest %>% select(-c(Study, Completeness, Contamination) ) %>% as.data.frame() %>% select( c("Sample" , Sign)) %>% column_to_rownames("Sample") %>% pheatmap::pheatmap(annotation_row = AnnotationColor, show_rownames = F)
  } else {
    #2. Only with 'Clade Low'
    UnirefTest %>% filter(Clade != "Other") %>% select(-c(Study, Completeness, Contamination) ) %>% as.data.frame() %>% select( c("Sample" , Sign)) %>% column_to_rownames("Sample") %>% t()  %>% pheatmap::pheatmap(annotation_col=AnnotationColor, show_colnames = F) 
  }
  return(list(Results_test,UnirefTest ))
}  


  
Enrichment_per_gene = Do_enrichment_uniref90(Feature = "Gene")
Enrichment_per_product = Do_enrichment_uniref90(Feature = "Product")


Enrichment_per_gene_all = Do_enrichment_uniref90(Feature = "Gene", Subset = F)
Enrichment_per_gene_all[[1]] %>% write_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/XuQ_analysis/Erichment_BMI_bromii.tsv")

Enrichment_per_gene_all[[1]] %>%  filter(Feature == "CladeOther") -> EnrichmentOther
Enrichment_per_gene_all[[1]] %>%  filter(Feature == "CladeLow2") -> EnrichmentLow

Enrichment_per_gene_all_fisher = Do_enrichment_uniref90(Feature = "Gene", Subset = F, Test = "fisher")
Enrichment_per_gene_all_fisher[[1]] %>% write_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/XuQ_analysis/Erichment_BMI_bromii_fisher.tsv")
Enrichment_per_gene_all_fisher[[1]] %>% left_join( UniRef90 %>% select(Gene, Product) %>% distinct() ) %>% arrange(P) -> Enrichment_per_gene_all_fisher[[1]]
ggplot(Enrichment_per_gene_all_fisher[[1]], aes(x=log10(Odds_ratio), y=-log10(P),col=FDR<0.05 ) ) + geom_point() + theme_bw()  + facet_wrap(~Comparison, scales= "free")

Enrichment_per_gene_all_fisher[[1]] %>% distinct(Comparison, Gene, .keep_all=T) %>% group_by(Comparison, FDR<0.05) %>% summarise(n())
Enrichment_per_gene_all_fisher[[1]] %>% distinct(Comparison, Gene, .keep_all=T) %>% filter(FDR<0.05) %>% group_by(Gene) %>% summarise(N = n()) %>% filter(N>1) #Diff enriched in both

#how many have annotation
#cazymes
Cazymes = read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Match_anpan/For_enrichment/Enrichment_tables/Cazymes.rds")
Enrichment_per_gene_all_fisher[[1]] %>% left_join(reaction_to_pathway) %>% group_by(is.na(Caz)) %>% summarise(n())
#EC
read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/For_enrichment/UniRef90/ec.rds") -> ECs
ECs_tibble <- enframe(ECs, name = "EC", value = "UniRef90")
ECs_tibble %>% separate_rows(UniRef90, sep="\t") -> ECs_tibble
Enrichment_tables = read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/For_enrichment/metacyc_pathways.rds")
Translate_EC_to_reaction = read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis//For_enrichment/EC_to_reaction.tsv.gz")
Enrichment_per_gene_all_fisher[[1]] %>% left_join(ECs_tibble %>% rename(Gene = UniRef90), relationship = "many-to-many" ) -> merged_ec
merged_ec %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #382/2090 annotated
#KEGG
read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/For_enrichment/UniRef90/KO.rds") %>% stack() %>% as_tibble() %>% rename(Gene = values, KO = ind) -> gene_to_ko
Enrichment_per_gene_all_fisher[[1]] %>% left_join(gene_to_ko ) -> merged_ko
merged_ko %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #0/2090 annotated
#Go
read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/For_enrichment/UniRef90/Go.rds") -> go_terms
go_terms %>% stack() %>% as_tibble() %>% rename(Gene = values, Go = ind) -> gene_to_go

#enrichment of GO
Enrichment_per_gene_all_fisher[[1]] %>% dplyr::select(-Product) %>% distinct() %>% dplyr::mutate(Odds_ratio =  log10(Odds_ratio) ) %>%  dplyr::rename(`z value` = Odds_ratio )  -> Input5 
Association_go = Run_enrichment(Input5 %>% filter(Comparison == "Other")  ,  Table=go_terms ) -> Go_other
Go_other %>% as_tibble() %>% filter(padj < 0.05  )
Association_go = Run_enrichment(Input5 %>% filter(Comparison != "Other")  ,  Table=go_terms ) -> Go_small
Go_small %>% as_tibble() %>% filter(padj < 0.05  )
#enrichment of pathways
Enrichment_per_gene_all_fisher[[1]]  %>% filter(Comparison == "Other") %>% mutate(p.adjust(P , "fdr") ) %>% rename(UniRef90 = Gene) %>% left_join(ECs_tibble) %>% rename(ec = EC) %>% left_join(Translate_EC_to_reaction) %>%  drop_na()  %>%  drop_na() -> Input5
Association_pathway = Run_enrichment(Input5 %>% rename(Gene = Reaction, `z value` = Odds_ratio ) ,  Table=Enrichment_tables )






Enrichment_per_gene_all_fisher[[1]] %>% left_join(gene_to_go ) -> merged_go
merged_go %>% drop_na() %>% group_by(Gene) %>% summarise(n()) #1,537/2090 annotated
Results_go = tibble()
for (Goo in unique(filter(merged_go, FDR<0.05)$Go)){
  if (is.na(Goo) | Goo == "NA" ){ next }
  merged_go %>% filter(Go == Goo) %>% group_by(FDR<0.05) %>% filter(Odds_ratio > 1) %>% summarise(N=n()) -> Term_up
  merged_go %>% filter(Go == Goo) %>% group_by(FDR<0.05) %>% filter(Odds_ratio < 1) %>% summarise(N=n()) -> Term_down
  
  merged_go %>% filter(Go != Goo) %>% group_by(FDR<0.05) %>% summarise(N=n()) -> Universe
  
  for (Term_n in c(1, 2)){
  Term = list(Term_up, Term_down)[[Term_n]]
  M = matrix(nrow = 2, ncol = 2)
  M[1,1] =  ifelse(T %in% Term$`FDR < 0.05`,filter(Term, `FDR < 0.05`==T)$N, 0)
  M[2,1] =  ifelse(F %in% Term$`FDR < 0.05`,filter(Term, `FDR < 0.05`==F)$N, 0) 
  M[1,2] =  ifelse(T %in% Universe$`FDR < 0.05`,filter(Universe, `FDR < 0.05`==T)$N, 0)
  M[2,2] =  ifelse(F %in% Universe$`FDR < 0.05`,filter(Universe, `FDR < 0.05`==F)$N, 0 )  
  fisher.test(M+1) -> Fisher_result
  tibble(Estimate = Fisher_result$estimate, P=Fisher_result$p.value, Direction=c("Up", "Down")[Term_n], Go_term = Goo  ) %>% rbind(Results_go, . ) -> Results_go
  }
} 
Results_go %>% filter(Direction == "Up") %>% mutate(FDR=p.adjust(P, "fdr")) %>% filter(FDR<0.05) -> Check_up
Results_go %>% filter(Direction == "Down") %>% mutate(FDR=p.adjust(P, "fdr")) %>% filter(FDR<0.05) -> Check_down
library(GO.db)
select(GO.db, keys = Check_up$Go_term , columns = c("GOID", "TERM", "ONTOLOGY", "DEFINITION"))
select(GO.db, keys = Check_down$Go_term , columns = c("GOID", "TERM", "ONTOLOGY", "DEFINITION"))

Explore_Go_enriched = function(Check_go){
  library(rrvgo)
  simMatrix <- calculateSimMatrix(Check_go$Go_term , orgdb="org.EcK12.eg.db", ont="MF", method="Rel")
  scores <- setNames(-log10(Check_go$P), Check_go$Go_term)
  reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.EcK12.eg.db")
  heatmapPlot(simMatrix, reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=6) -> Plot1
  scatterPlot(simMatrix, reducedTerms) ->Plot2
  treemapPlot(reducedTerms) -> Plot3
  return(list(Plot1, Plot2, Plot3))
}


##Heatmap
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  tiff(filename, width=width, height=height, units= "in", res=300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

left_join(Enrichment_per_gene, UniRef90 %>% select(Gene, Product) ) -> Enrichment_per_gene
Enrichment_per_gene %>% mutate(FDR = p.adjust( `Pr(>|z|)`, "fdr") ) %>% filter(FDR<0.05) %>% group_by(Product) %>% summarise(n())

Enrichment_per_gene_all[[1]] %>% mutate(FDR= p.adjust( `Pr(>|z|)`, "fdr") ) %>% filter(FDR<0.05) -> Sign
Enrichment_per_product[[2]]  %>% select(-c(Study, Completeness, Contamination) ) %>% as.data.frame() %>% select( c("Sample" , Sign$Gene)) %>% column_to_rownames("Sample") %>% t()  %>% pheatmap::pheatmap(annotation_col=AnnotationColor, show_colnames = F, fontsize_row = 9) 

Enrichment_per_gene_all_fisher[[1]] %>% filter(FDR<0.005) -> Sign
Enrichment_per_gene_all_fisher[[2]] %>% select(Sample, Clade ) %>% as.data.frame()   %>% column_to_rownames("Sample") %>% drop_na() -> AnnotationColor
Enrichment_per_gene_all_fisher[[2]] %>% drop_na() %>% select(-c(Study, Completeness, Contamination) ) %>% as.data.frame() %>% select( c("Sample" , Sign$Gene)) %>% column_to_rownames("Sample") %>% t()  %>% pheatmap::pheatmap(annotation_col=AnnotationColor, show_colnames = F, show_rownames = F) 

Enrichment_per_gene_all_fisher[[2]] %>% ggplot(aes(x=Completeness, fill=Clade )) + geom_density(alpha=0.4) + theme_bw()
Enrichment_per_gene_all_fisher[[2]] %>% ggplot(aes(x=Contamination, fill=Clade )) + geom_density(alpha=0.4) + theme_bw()


###Upset comparing all results
Enrichment_per_gene_all_fisher[[1]] %>% filter(Comparison=="Other") %>% mutate(FDR2 = p.adjust(P, "fdr") ) %>% filter(FDR2<0.05) %>% arrange(P) %>% mutate(Type = "Fisher_other", Presence = 1) %>% select(Type, Gene, Presence) %>% spread(Gene, Presence) -> Fisher_other
Enrichment_per_gene_all_fisher[[1]] %>% filter(Comparison!="Other") %>% mutate(FDR2 = p.adjust(P, "fdr") ) %>% filter(FDR2<0.05) %>% arrange(P)  %>% mutate(Type = "Fisher_low", Presence = 1) %>% select(Type, Gene, Presence) %>% spread(Gene, Presence)  -> Fisher_low
Enrichment_per_gene_all[[1]] %>% filter(Feature == "CladeLow2") %>% mutate(FDR2 = p.adjust(`Pr(>|z|)`, "fdr") ) %>% filter(FDR2<0.05) %>% arrange(`Pr(>|z|)`)  %>% mutate(Type = "GLM_other", Presence = 1) %>% select(Type, Gene, Presence) %>% spread(Gene, Presence)  -> GLM_other
Enrichment_per_gene_all[[1]] %>% filter(Feature != "CladeLow2") %>% mutate(FDR2 = p.adjust(`Pr(>|z|)`, "fdr") ) %>% filter(FDR2<0.05) %>% arrange(`Pr(>|z|)`)  %>% mutate(Type = "GLM_low", Presence = 1) %>% select(Type, Gene, Presence) %>% spread(Gene, Presence)  -> GLM_low

full_join(Fisher_other, Fisher_low) %>% full_join(GLM_other) %>% full_join(GLM_low) %>% as.data.frame() %>% column_to_rownames("Type") -> Sets_Compare
Sets_Compare[is.na(Sets_Compare)] = 0
UpSetR::upset(Sets_Compare %>% t()  %>% as.data.frame()) #, sets.x.label = "SGB tree-Phenotype\nsupported assocations", mainbar.y.label = "Number intersected supported associatons", text.scale = 3, nsets = 6, nintersects = 15 ) -> Plot2





#######ENRICHMENT OF ASSOCIATIONS##############
Run_enrichment = function(Summary_stats, Table,Plot=F){
  Summary_stats %>%  arrange(desc(`z value`))-> R1
  #Prepare ranks
  ranks = R1$`z value`
  names(ranks) = as.character(R1$Gene)
  #Run GSEA
  fgseaRes <- fgsea::fgsea(Table, ranks)
  
  if (Plot == T){
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    fgsea::plotGseaTable(Table[topPathways], ranks, fgseaRes,  gseaParam=0.5)
  }
  
  return(fgseaRes)
}
Create_list_function = function(File){
  #Opens non-header file in whcih the first column should be the name of a list item, and there is a subsequent number of unidetermined column numbers that should become the vector of item for each key
  dat <- readLines(File)
  dat <- strsplit(dat, "\t")
  
  lapply(dat, function(entry) {
    entry[1] # Remove the first value
  }) %>% unlist() -> Names
  
  data_list <- lapply(dat, function(entry) {
    entry[-1]  # Remove the first value
  })
  names(data_list) = Names
  
  return(data_list)
}

#Is there a cazyme enrichment?
Cazymes = read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Match_anpan/For_enrichment/Enrichment_tables/Cazymes.rds")
reaction_to_pathway <- stack(Cazymes)
as_tibble(reaction_to_pathway) %>% rename(Gene=values, Caz=ind) -> r_to_p

Make_analysis = function(SummaryStats, reaction_to_pathway = r_to_p ){
  
  #Match
  SummaryStats %>% left_join(reaction_to_pathway) -> SummaryStats
  
  #Get for enrichment
  SummaryStats %>% filter( is.na(Caz) ) -> Not_caz
  SummaryStats %>% filter(! is.na(Caz) ) -> caz
  if (dim(Not_caz)[1] ==0 | dim(caz)[1] == 0){ return() }
  Cazymes_enrichment = list( "Not_cazymes" = Not_caz$Gene, "Cazymes" = caz$Gene)
  
  Run_enrichment(SummaryStats, Cazymes_enrichment) -> Result
  return(Result)
}

Caz= Make_analysis(Enrichment_per_gene)
Caz2 = Make_analysis(Enrichment_per_gene_all[[1]] %>% filter(Feature == "CladeOther") )
Caz %>% filter(!pathway == "Not_cazymes") %>% select(-c(leadingEdge, File, padj)) %>% arrange(pval) %>% mutate(FDR = p.adjust(pval, "fdr")) -> Caz

#Is there a pathway enrichment? --> Need enrichment of EC

#Make UniRef90 into EC
read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/For_enrichment/UniRef90/ec.rds") -> ECs
ECs_tibble <- enframe(ECs, name = "EC", value = "UniRef90")
ECs_tibble %>% separate_rows(UniRef90, sep="\t") %>% filter(UniRef90 %in% Enrichment_per_gene$Gene ) -> ECs_tibble


Enrichment_tables = read_rds("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/For_enrichment/metacyc_pathways.rds")
Translate_EC_to_reaction = read_tsv("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis//For_enrichment/EC_to_reaction.tsv.gz")
#1. Make ECs to reactions
Enrichment_per_gene %>% rename(UniRef90 = Gene) %>% left_join(ECs_tibble) %>% rename(ec = EC) %>% left_join(Translate_EC_to_reaction) %>%  drop_na() -> Input #only check ECs where we have annotation of the reaction and pathway
Enrichment_per_gene_all[[1]]  %>% filter(Feature == "CladeOther") %>% mutate(p.adjust(`Pr(>|z|)` , "fdr") ) %>% rename(UniRef90 = Gene) %>% left_join(ECs_tibble) %>% rename(ec = EC) %>% left_join(Translate_EC_to_reaction) %>%  drop_na()  %>%  drop_na() -> Input2
Enrichment_per_gene_all[[1]] %>% filter(Feature != "CladeOther") %>% rename(UniRef90 = Gene) %>% left_join(ECs_tibble) %>% rename(ec = EC) %>% left_join(Translate_EC_to_reaction) %>%  drop_na()  %>%  drop_na() -> Input3
Enrichment_per_gene_all_fisher[[1]]  %>% filter(Comparison == "Other") %>% mutate(p.adjust(P , "fdr") ) %>% rename(UniRef90 = Gene) %>% left_join(ECs_tibble) %>% rename(ec = EC) %>% left_join(Translate_EC_to_reaction) %>%  drop_na()  %>%  drop_na() -> Input4
Enrichment_per_gene_all_fisher[[1]]  %>% filter(Comparison != "Other") %>% mutate(p.adjust(P , "fdr") ) %>% rename(UniRef90 = Gene) %>% left_join(ECs_tibble) %>% rename(ec = EC) %>% left_join(Translate_EC_to_reaction) %>%  drop_na()  %>%  drop_na() -> Input5


#2. Look for pathway enrichment
Association = Run_enrichment(Input %>% rename(Gene = Reaction) ,  Table=Enrichment_tables )
Association %>% as_tibble() %>% arrange(pval)


Association2 = Run_enrichment(Input2 %>% rename(Gene = Reaction) ,  Table=Enrichment_tables )
Association2 %>% as_tibble() %>% filter(padj < 0.05)

r1 = Input2 %>%  arrange(desc(`z value`))
ranks = r1$`z value` ; names(ranks) = as.character(Input2$Reaction)

plotEnrichment(Enrichment_tables[["PWY-6897"]], ranks) + labs(title="PWY-6897") + theme( text = element_text(size = 21),  plot.title = element_text(size = 26), axis.title.x = element_text(size = 21),  axis.title.y = element_text(size = 21)  )
plotEnrichment(Enrichment_tables[["PWY-7400"]], ranks) + labs(title="PWY-7400") + theme( text = element_text(size = 21),  plot.title = element_text(size = 26), axis.title.x = element_text(size = 21),  axis.title.y = element_text(size = 21)  )

fora(Enrichment_tables, filter(Input2, FDR<0.05 & Estimate>0)$ec , unique(Input2$ec) , minSize = 1, maxSize = Inf)
fora(Enrichment_tables, filter(Input2, FDR<0.05 & Estimate<0)$ec , unique(Input2$ec) , minSize = 1, maxSize = Inf)

fora(Enrichment_tables, filter(mutate(Input, FDR=p.adjust(`Pr(>|z|)`,"fdr")), FDR<0.05 & Estimate>0)$ec , unique(Input2$ec) , minSize = 1, maxSize = Inf)
fora(Enrichment_tables, filter(mutate(Input, FDR=p.adjust(`Pr(>|z|)`,"fdr")), FDR<0.05 & Estimate<0)$ec , unique(Input2$ec) , minSize = 1, maxSize = Inf)


fora(Enrichment_tables, filter(Input4,, FDR<0.05 & Odds_ratio>1)$ec , unique(Input4$ec) , minSize = 1, maxSize = Inf)
Association3 = Run_enrichment(Input4 %>% rename(Gene = Reaction, `z value` = Odds_ratio ) ,  Table=Enrichment_tables )

plotEnrichment(Enrichment_tables[["PWY-6897"]], ranks) + labs(title="PWY-6897") + theme( text = element_text(size = 21),  plot.title = element_text(size = 26), axis.title.x = element_text(size = 21),  axis.title.y = element_text(size = 21)  )
plotEnrichment(Enrichment_tables[["PWY-7400"]], ranks) + labs(title="PWY-7400") + theme( text = element_text(size = 21),  plot.title = element_text(size = 26), axis.title.x = element_text(size = 21),  axis.title.y = element_text(size = 21)  )


#Do the same in other clade
Clade_low_vs_high = function(SGB="t__SGB4933_group", Tree_path="/mnt/project/Symlink_phylo/IQtree.t__SGB4933_group.TreeShrink.tre" ){
  Uniref = paste0("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Gene_presence/",SGB, ".tsv")
  #Get stats from anpan. Need to save them first.
  Anpan_stats = paste0("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/XuQ_analysis/BMI_",SGB, ".tsv") -> SGB_BMI
  #Read tree
  Tree = Process_tree(Tree_path)
  #read
  read_tsv(Uniref) -> UniRef90
  read_tsv(SGB_BMI) -> SGB_BMI
  
  
  #Filters
  UniRef90 %>% filter(Sample %in% SGB_BMI$sample_id ) -> UniRef90
  SGB_BMI %>% mutate(MAG_available = ifelse(sample_id %in% UniRef90$Sample, T, F) ) -> SGB_BMI
  keep.tip(Tree,SGB_BMI$sample_id) -> Tree
  
  #annotation
  SGB_BMI %>% mutate(BMI_group = ifelse(BMI<18.5, "BMI<18.5", ifelse(BMI>=18.5 & BMI<25, "18.5≤BMI<25", ifelse(BMI>=25 & BMI<30, "25≤BMI<30", "30≤BMI" ) ) ) ) -> SGB_BMI
  SGB_BMI %>% mutate(Age_group = ifelse(Age ==0, "<1", ifelse(Age < 15, "(0-15)", ifelse( Age < 30, "[15, 30)", ifelse( Age < 70, "[30, 70)", ">=70") )))) %>% mutate(Age_group = factor(Age_group, levels=c( "<1", "(0-15)", "[15, 30)", "[30, 70)", ">=70"  ))) -> SGB_BMI
  
  #Take subset
  set.seed(12)
  sample_n = ifelse(dim(SGB_BMI)[1] > 1000, 1000, dim(SGB_BMI)[1])
  sample(SGB_BMI$sample_id, sample_n) -> SubSet
  SGB_BMI %>% filter(sample_id %in% SubSet) -> SGB_BMI2
  keep.tip(Tree, SubSet) -> Tree2
  
  #Find ancestor clades "high" and "low" from subset
  SGB_BMI2 %>% arrange(phylo_effect_median) %>% head(n=100)  %>% select(sample_id) %>% as_vector() ->   leafs_ids   
  leaf_nodes <- match(leafs_ids, Tree2$tip.label)
  leaf_nodes = seq(387, 403)
  #leaf_nodes = seq(245, 263) $t__SGB4910
  #leaf_nodes = seq(87,123) #for SGB15254/Oscillibacter_sp_ER4
  #leaf_nodes = seq(67, 75) #for t__SGB15318_group
  common_ancestor1 <- getMRCA(Tree2,leaf_nodes) #
  Names1 = Tree2$tip.label[ leaf_nodes]
  
  SGB_BMI2 %>% arrange(desc(phylo_effect_median)) %>% head(n=20)  %>% select(sample_id) %>% as_vector() ->   leafs_ids   
  leaf_nodes <- match(leafs_ids, Tree2$tip.label)
  leaf_nodes = leaf_nodes = seq(260,265)
  #leaf_nodes = seq(566, 577) #t__SGB4910
  #leaf_nodes = seq(530, 555) #for SGB15254/Oscillibacter_sp_ER4
  #leaf_nodes = seq(902,914) #for t__SGB15318_group
  common_ancestor2 <- getMRCA(Tree2, leaf_nodes )
  Names2 = Tree2$tip.label[leaf_nodes]
  
  
  #Make plot
  
  ggtree(Tree2, layout = "fan",  open.angle=15, size=0.1 ) %<+% SGB_BMI2 -> p
  p + geom_fruit( geom="geom_tile", mapping = aes(fill=BMI_group), width=0.03,offset=0.1  ) + scale_fill_manual(values = c("BMI<18.5" = "blue",  "18.5≤BMI<25"= "grey", "25≤BMI<30" = "#ff7f0e", "30≤BMI"= "red"), na.value="white" ) +   labs(fill = "BMI") + 
    new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Age_group), width=0.03,offset=0.1  ) + scale_fill_manual(values = c("<1" = "#5E1DB5",  "(0-15)"= "#5A3C82", "[15, 30)" = "#0E27E8", "[30, 70)"= "#EBB649",  ">=70"="#B5731D" ), na.value="white" ) +   labs(fill = "Age") +
    new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
    geom_balance(node=common_ancestor1  , fill="skyblue", color=NA, alpha=0.3) + geom_balance(node=common_ancestor2  , fill="#FF7F7F", color=NA, alpha=0.3)
  
  
  #For testing
  
  
  #Find ancestor clades "high" and "low"
  common_ancestor3 <- getMRCA(Tree,  match(Names1, Tree$tip.label) )
  common_ancestor4 <- getMRCA(Tree,  match(Names2, Tree$tip.label) )
  
  
  Tree$tip.label[getDescendants(Tree, node = common_ancestor3)] -> Low_BMI
  Tree$tip.label[getDescendants(Tree, node = common_ancestor4)] -> High_BMI
  
  SGB_BMI %>% mutate( Clade = ifelse(sample_id %in% Low_BMI, "Low", ifelse(sample_id %in% High_BMI, "High", NA)  )   ) %>% drop_na() -> ToTest
  Phenotypes %>% select(ID_anal, study_name, Sex) %>% rename(sample_id = ID_anal ) %>% left_join(ToTest, . ) -> ToTest
  
  
  ToTest %>% filter(Country %in% c("USA")) %>%  ggplot(aes(x=Clade, y=BMI)) + geom_boxplot(outlier.shape = NA, aes(fill=Clade)) + theme_bw() + geom_sina(alpha=0.5 ) + facet_wrap(~Country , scales="free") +  
    scale_fill_manual(values = c("High"= "#FF7F7F", "Low" = "skyblue") ) + theme(text = element_text(size = 14), axis.text.y = element_blank()) + coord_flip()
  ToTest %>% filter(Country %in% c( "USA") ) %>%  ggplot(aes(x=Clade, y=BMI, col=study_name  )) + geom_boxplot(outlier.shape = NA) + theme_bw() + geom_sina(alpha=0.5 ) + facet_wrap(~Country , scales="free") +  
    theme(text = element_text(size = 14) ) + coord_flip() + scale_color_manual(values =c25)
  
  
  
  ToTest %>% filter(Country %in% c("NLD")) %>% drop_na() %>% lmerTest::lmer( BMI ~ Clade + Age + Sex + (1|Country) + (1|Country:study_name), . ) %>% summary()
  
  
  #Check for MAGs
  ToTest %>% group_by(Clade, MAG_available) %>% summarise(n()) %>% print()

  #If enough ToTest %>% filter(MAG_available == T)
  ToTest %>% filter(MAG_available == T)
  UniRef90 %>% filter(Sample %in% filter(ToTest, MAG_available == T)$sample_id ) -> UniRef_test
  UniRef_test %>% select(Sample, Product) %>% mutate(N=1) %>%  spread(Product, N)
  
  
  UniRef_test %>% distinct(Gene, Sample) %>% group_by(Gene) %>% summarise(N = n()) %>% arrange(desc(N)) %>% filter(! (N > 9 | N < 3 ) ) -> Keep
  UniRef_test %>% filter(Gene %in% Keep$Product) %>% distinct(Sample, Gene) %>% mutate(N=1) %>% spread(Gene, N) -> Wide_df
  Wide_df[is.na(Wide_df)] = 0
  ToTest %>% filter(MAG_available == T) %>% as.data.frame() %>% column_to_rownames("sample_id") %>% select(Clade) -> Annotation
  Wide_df %>% as.data.frame() %>% column_to_rownames("Sample") %>% pheatmap::pheatmap(. , annotation_row = Annotation, annotation_colors = list(Group = c(Low = "red", High = "blue")) )
  
  
  Uniref = paste0("/mnt/project/Make_Associations/Functional_enrichment/MAG_analysis/Gene_presence/Cazymes/",SGB, ".tsv")
  read_tsv(Uniref) -> CAZ
  
  
  
}


SGB = "t__SGB15318_group"


#########################
#Tree comparison########
########################

#Human tree#
read.tree("/mnt/project/Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Data/fst_min_evolution.nwk") -> Human_tree
read.tree("/mnt/project/Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Data/fst_upgma.nwk") -> Human_tree


Compare_trees = function(Human_tree, bacteria_file ){

  #Bacterial tree#
  read_tsv(bacteria_file) -> Bug_tree
  Bug_tree %>%  separate(Country_ID, into = c("Country1", "Country2"), sep = "-") %>% mutate(Median_distance_by_country = ifelse(Country1==Country2, 0 , Median_distance_by_country) ) %>% 
    filter( (Country1 %in% Human_tree$tip.label) & (Country2 %in% Human_tree$tip.label) ) -> long_distance
  other_side = tibble()
  for (Entry in seq(dim(long_distance)[1]) ){
    long_distance[Entry,] -> Entry
    if (Entry$Country2 == Entry$Country1){ next }
    other_side %>% rbind( tibble(Country1 = Entry$Country2, Country2=Entry$Country1, Median_distance_by_country= Entry$Median_distance_by_country ) )  -> other_side
    
  }
  long_distance %>% rbind(other_side) %>%  pivot_wider(names_from = Country2, values_from = Median_distance_by_country) -> Wide_distance
  Wide_distance[is.na(Wide_distance)] = 0
  
  
  # Remove the Country1 column to get only the distance values
  distance_matrix <- as.matrix(Wide_distance[, -1])
  # Set the row and column names of the matrix
  rownames(distance_matrix) <- Wide_distance$Country1
  colnames(distance_matrix) <-  colnames(Wide_distance[,-1])
  distance_matrix %>% as.dist() %>%  ape::njs() -> Bug_tree2
  ggtree(Bug_tree2) + geom_tiplab()
  
  
  keep.tip(Human_tree, unique(c(long_distance$Country1,long_distance$Country2))) -> Human_tree2
  distance_human = cophenetic.phylo( Human_tree2 ) 
  
  distance_matrix <- distance_matrix[rownames(distance_human), colnames(distance_human)]
  
  
  mantel.test(distance_matrix, distance_human) -> Testresults
  
  return( list(Testresults, Bug_tree2, Human_tree2  ))

}

Patt = "_country.tsv"
matching_files <- list.files("/mnt/project/Strain_sharedness/Distance_countryNcontinent/", pattern = paste0("*", Patt) , full.names = TRUE)
Summary_stats_tree2tree = tibble()
for (File in matching_files){
  str_split(File, "/")[[1]] -> SGB ;  SGB = SGB[length(SGB)] ; str_replace(SGB, Patt, "") -> SGB
  Compare_trees(Human_tree, File) -> result
  Summary_stats_tree2tree = rbind(Summary_stats_tree2tree, tibble(SGB= SGB, rho = result[[1]]$z.stat, P= result[[1]]$p ) )
  
}

Summary_stats_tree2tree %>% arrange(P) %>% mutate(FDR = p.adjust(P, "fdr") ) %>% filter(FDR<0.05) %>% arrange(desc(rho))

Compare_trees(Human_tree, matching_files[grepl("t__SGB6750" , matching_files)] ) -> check


check[[2]]$edge.length = 1

ggtree( check[[2]],  ) + geom_tiplab()
ggtree( check[[3]],  ) + geom_tiplab()

comparePhylo(check[[2]], check[[3]], plot = T, force.rooted = FALSE, use.edge.length = FALSE)



T1 <- ggtree(check[[2]]) +    
  theme_tree2(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) +   
  geom_tiplab()   
T2 <- ggtree(check[[3]]) +   
  theme_tree2(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) +   
  geom_tiplab(hjust =1) +   
  scale_x_reverse()     
d1 = T1$data[T1$data$isTip,]  
d1$x[] = 1  
d2 = T2$data[T2$data$isTip,]  
d2$x[] = 2  

TTcon <- rbind(d1, d2)  

T1 = ggtree(check[[2]] ) +    
  theme_tree2(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) +   
  geom_tiplab(align = T) +  xlim(0,2.5)  
T2 = ggtree(check[[3]]) +   
  theme_tree2(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) +   
  geom_tiplab(hjust =1, align = T) +   
  scale_x_reverse( limits = c(0.035, 0))   

L1 = ggplot(TTcon, aes(x = x, y = y, colour = label, group = label)) + geom_line() +   
  theme_void() + theme(legend.position="none", plot.margin = unit(c(1,0,1,0),"cm"))  

cowplot::plot_grid(T1, L1 ,T2, nrow = 1, align = "hv")



################################################
############IBD plot###########################
###############################################

#1. Get all  IBD-related scores
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Europe.csv") %>% filter(Phenotype == "IBD") -> All_Europe
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Asia.csv") %>% filter(Phenotype == "IBD")  -> All_Asia
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica.csv") %>% filter(Phenotype == "IBD") -> All_NorthAmerica

SGBs_for_tree = unique(c(All_Europe$SGB, All_Asia$SGB, All_NorthAmerica$SGB ))


read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Europe_filtered.csv") %>% filter(Phenotype == "IBD") -> sign_Europe
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Asia_filtered.csv") %>% filter(Phenotype == "IBD") -> sign_Asia
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica_filtered.csv") %>% filter(Phenotype == "IBD") -> sign_NorthAmerica

All_Europe %>% mutate(Supported = ifelse(SGB %in% sign_Europe$SGB, T, F) ) %>% select(SGB, elpd_diff, Supported)  -> All_Europe
All_Asia %>% mutate(Supported = ifelse(SGB %in% sign_Asia$SGB, T, F) ) %>% select(SGB,elpd_diff, Supported) -> All_Asia
All_NorthAmerica %>% mutate(Supported = ifelse(SGB %in% sign_NorthAmerica$SGB, T, F) ) %>% select(SGB,elpd_diff, Supported) -> All_NorthAmerica

##Associations when controlling for prevalence
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_Europe.csv")%>% filter(Phenotype == "IBD") -> prevalence_Europe
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_Asia.csv")%>% filter(Phenotype == "IBD") -> prevalence_Asia
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_NorthAmerica.csv")%>% filter(Phenotype == "IBD") -> prevalence_NorthAmerica

read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_Europe_filtered.csv")%>% filter(Phenotype == "IBD") -> sign_Europe
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_Asia_filtered.csv")%>% filter(Phenotype == "IBD") -> sign_Asia
read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_NorthAmerica_filtered.csv")%>% filter(Phenotype == "IBD") -> sign_NorthAmerica

prevalence_Europe %>% mutate(Supported = ifelse(SGB %in% sign_Europe$SGB, T, F) ) %>% select(SGB,elpd_diff, Supported)  -> prevalence_Europe
prevalence_Asia %>% mutate(Supported = ifelse(SGB %in% sign_Asia$SGB, T, F) ) %>% select(SGB,elpd_diff, Supported) -> prevalence_Asia
prevalence_NorthAmerica %>% mutate(Supported = ifelse(SGB %in% sign_NorthAmerica$SGB, T, F) ) %>% select(SGB,elpd_diff, Supported)  -> prevalence_NorthAmerica

##Association in Dutch
read_csv("/mnt/project/Make_Associations/Association/Results_DutchIBD.csv") %>%  filter(Phenotype == "IBD")   -> Dutch
read_csv("/mnt/project/Make_Associations/Association/Results_DutchIBD_filtered.csv") %>% filter(Phenotype == "IBD")  -> sign_Dutch
Dutch %>% mutate(Supported = ifelse(SGB %in% sign_Dutch$SGB, T, F) ) %>% select(SGB,elpd_diff, Supported)  -> Dutch

##Association CD/UC
read_csv("/mnt/project/Make_Associations/Association/Results_UC_CD.csv") %>% mutate(ID = paste0(Phenotype, SGB) )    -> UCCD
read_csv("/mnt/project/Make_Associations/Association/Results_UC_CD_filtered.csv") %>% mutate(ID = paste0(Phenotype, SGB) )   -> sign_UCCD
UCCD %>% mutate(Supported = ifelse(ID %in% sign_UCCD$ID, T, F) )   %>% select(SGB,elpd_diff, Supported, Phenotype)  -> UCCD



#2. get SGB tree
str_replace(SGBs_for_tree, "t__", "") -> SGBs_for_tree
Info %>% filter(SGB %in% SGBs_for_tree) -> SGBs_for_tree

read.tree("/mnt/project/Make_Associations/Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk") -> SGB_tree_all
#reroot
SGB_tree_all %>% root(., outgroup = 2066) -> SGB_tree_all
keep.tip(SGB_tree_all, SGBs_for_tree$SGB_name) -> SGB_tree_IBD

SGBs_for_tree %>% mutate(SGB = paste0("t__", SGBs_for_tree$SGB) ) %>% select(SGB_name, SGB, Taxonomy, Phylum ) -> ForPlot #%>% left_join(. ,  All_Europe  ) %>% left_join(. , All_Asia, by="SGB",  suffix=c("_Europe", "_Asia") ) %>%  left_join(. , All_NorthAmerica, by="SGB"  ) -> ForPlot
  
ForPlot %>% mutate(Genus = str_split(Taxonomy, "\\|", simplify = TRUE)[, 6]) %>% mutate(Genus = str_split(Genus, "__", simplify = TRUE)[, 2]) -> ForPlot
ForPlot %>% mutate(ID_name = paste0( str_replace(SGB, "t__", ""),"\n",str_replace(Genus, "_unclassified", "") ) )  -> ForPlot

  
#3. Make plot

#With barplots
rbind(rbind(All_Europe %>% mutate(Continent="Europe"), All_Asia%>% mutate(Continent="Asia")), All_NorthAmerica%>% mutate(Continent="North_America")) -> ContinentSratifiedInfo
rbind(rbind(prevalence_Europe %>% mutate(Continent="Europe"), prevalence_Asia%>% mutate(Continent="Asia")), prevalence_NorthAmerica%>% mutate(Continent="North_America")) -> ContinentSratified_prevalence_Info
left_join(tibble(SGB_name = SGB_tree_IBD$tip.label), ForPlot)$SGB -> SGB_tree_IBD$tip.label

ggtree(SGB_tree_IBD,  layout = "fan",  open.angle=15, size=0.1) +  geom_fruit(data=left_join(ContinentSratifiedInfo, ForPlot) ,  geom=geom_bar, mapping=aes(x=elpd_diff, y=SGB, fill=Continent, color=Supported), stat="identity", orientation="y", offset=0.2) +  scale_color_manual(values = c("TRUE"="black", "FALSE"="grey" ) ) +
  geom_fruit(data=left_join(ContinentSratified_prevalence_Info, ForPlot) ,  geom=geom_bar, mapping=aes(x=elpd_diff, y=SGB, fill=Continent,color=Supported), stat="identity", orientation="y") +
  scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) +
  new_scale_fill() +
  geom_fruit(data=left_join(Dutch, ForPlot) ,  geom=geom_bar, mapping=aes(x=elpd_diff, y=SGB, color=Supported), stat="identity", orientation="y", , fill="lightblue") +
  geom_fruit(data=left_join(UCCD, ForPlot) ,  geom=geom_bar, mapping=aes(x=elpd_diff, y=SGB, fill=Phenotype, color=Supported ), stat="identity", orientation="y" ) + scale_fill_manual(values=c("CrohnsDisease" = "brown", "UlcerativeColitis" = "darkturquoise")) -> p
p %<+%  select(ForPlot, -SGB_name) +
   new_scale_color() +
  geom_tiplab(aes(label=ID_name, col=Phylum ) ,size=1.5) + scale_color_manual(values= c("Bacteroidetes"="green4","Firmicutes"="#E31A1C", "Actinobacteria"="dodgerblue2" ) ) -> p
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/IBD_associations_summary.pdf", p)  


#with tiles
ggtree(SGB_tree_IBD,  layout = "fan",  open.angle=15, size=0.1) +
  geom_fruit(data=left_join(All_Europe, ForPlot) ,  geom=geom_tile, mapping=aes(y=SGB, fill=elpd_diff, col=Supported), width=0.03,offset=0.2, size=1) + scale_color_manual(values = c("TRUE"="black", "FALSE"="white" ) )  +
  geom_fruit(data=left_join(All_Asia, ForPlot) ,  geom=geom_tile, mapping=aes(y=SGB, fill=elpd_diff, col=Supported), width=0.03,offset=0.2, size=1)   +
  geom_fruit(data=left_join(All_NorthAmerica, ForPlot) ,  geom=geom_tile, mapping=aes(y=SGB, fill=elpd_diff, col=Supported), width=0.03,offset=0.2, size=1)   +
  new_scale_fill() +
  geom_fruit(data=left_join(prevalence_Europe, ForPlot) ,  geom=geom_tile, mapping=aes(y=SGB, fill=elpd_diff, col=Supported), width=0.03,offset=0.2, size=1)   +
  geom_fruit(data=left_join(prevalence_Asia, ForPlot) ,  geom=geom_tile, mapping=aes(y=SGB, fill=elpd_diff, col=Supported), width=0.03,offset=0.2, size=1) +
  geom_fruit(data=left_join(prevalence_NorthAmerica, ForPlot) ,  geom=geom_tile, mapping=aes(y=SGB, fill=elpd_diff, col=Supported), width=0.03,offset=0.2, size=1)  + scale_fill_gradient(high = "lightblue", low = "darkblue") +
  new_scale_fill() +
  geom_fruit(data=left_join(Dutch, ForPlot) ,  geom=geom_tile, mapping=aes(y=SGB, fill=elpd_diff, col=Supported), width=0.03,offset=0.2, size=1)  + scale_fill_gradient(high = "#ff0000", low = "#300000") +
  new_scale_fill() +
  geom_fruit(data=left_join(UCCD %>% filter(Phenotype == "CrohnsDisease"), ForPlot),  geom=geom_tile, mapping=aes(y=SGB, fill=elpd_diff, col=Supported), width=0.03,offset=0.2, size=1) +
  geom_fruit(data=left_join(UCCD %>% filter(Phenotype == "UlcerativeColitis"), ForPlot) ,  geom=geom_tile, mapping=aes(y=SGB, fill=elpd_diff, col=Supported), width=0.03,offset=0.2, size=1) + scale_fill_gradient(high = "#0eff00", low = "#063b00") +
  new_scale_fill() + new_scale_color() -> p
  p %<+%  select(ForPlot, -SGB_name) +
  geom_tiplab(aes(label=ID_name, col=Phylum ) ,size=1) + scale_fill_manual(values= c("Bacteroidetes"="green4","Firmicutes"="#E31A1C", "Actinobacteria"="dodgerblue2" ) ) -> p
  

ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/IBD_associations_summary.pdf", p)  
  
  
  
######################
#####FST##############
######################
library("SNPRelate")
library(tidyverse)
Location = "/mnt/project/"
args <- commandArgs(trailingOnly = TRUE)
SGB = args[1]
VCF = paste0(Location,  "Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Generate_VCF/VCF/", SGB, ".vcf" )
GDS = paste0(Location,  "Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/GDS/", SGB, ".gds" )
Phenotype = paste0(Location,"Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv")
#Outputs
Bug_fst_file = paste0(Location, "Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Fst_matrices/",SGB,".csv")
Location_human = paste0(Location, "Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Fst_matrices/Human_fst.mat")
Association_summary = paste0(Location, "Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Association_Fst/", SGB, ".tsv")

#Generate GDS file. Method is the only one that allows non-biallelic entries
snpgdsVCF2GDS(VCF, GDS, method="copy.num.of.ref")
genofile <- snpgdsOpen(GDS)
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
#Get annotation of each sample
pop <- read_tsv(Phenotype) %>% dplyr::select(ID_anal, Country) %>% drop_na()
pop %>% filter(ID_anal %in% samples) %>% distinct(ID_anal, .keep_all=T) %>% as.data.frame() %>% column_to_rownames("ID_anal") -> pop
# Prepare to compute Fst
pop_ordered <- pop[samples,]
countries = unique(pop_ordered)
countries = countries[! is.na(countries) ]
fst_matrix <- matrix(nrow = length(countries), ncol = length(countries))
for (i in 1:length(countries)){
  for (j in 1:length(countries)){
    if (i == j) {
      fst_matrix[i,j] <- 0
    } else {
      c1 <- countries[i]
      c2 <- countries[j]
      if (is.na(c1) | is.na(c2) ){ next }
      
      # choose samples from the 2 countries used in Fst calculation
      flag <- pop_ordered %in% c(c1,c2) 
      samp.sel <- samples[flag]
      pop.sel <- pop_ordered[flag]
      table(pop.sel) -> Table_c
      if (sum(Table_c < 10 ) > 0 ) { fst_matrix[i,j] <- NA ; next }
      fst_res <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel), method = "W&H02" )
      fst_matrix[i,j] <- fst_res$Fst
    }
  }
}
row.names(fst_matrix) <- countries
colnames(fst_matrix) <- countries
snpgdsClose(genofile)
fst_matrix %>% apply(2, function(x){ x[!x==0] -> y ; sum(is.na(y))/length(y) }  ) -> Proportion_na
colnames(fst_matrix)[Proportion_na!=1] -> K
fst_matrix[K, K] -> fst_matrix

write.csv(fst_matrix, file = Bug_fst_file, row.names = T)

#Compare with human Fst
Human_matrix = read_tsv(Location_human) %>% as.data.frame() %>% column_to_rownames("...1")
intersect(rownames(fst_matrix), rownames(Human_matrix) ) -> countries_check
Human_matrix[countries_check,countries_check] -> Human_matrix
fst_matrix[countries_check,countries_check] -> Bug_matrix
Per = 10000
mantel( Human_matrix  , Bug_matrix,permutations = Per, method="spearman") -> Result_mantel

tibble(Rho=Result_mantel$z.stat,  P = Result_mantel$p, Permutations = 3000 ) %>% write_tsv(Association_summary)

#################################
#####Blongum tree with refs#####
###############################
library(tidyverse)
library(ape)
library(ggtree)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/Blongum_infant.tsv") -> longum_infant
#Process_tree("/mnt/project/Phylogenies/t__SGB17248_withLongumSubs/treeshrink/RAxML_bestTree.t__SGB17248.TreeShrink.tre") -> Tree
Process_tree("/mnt/project/Phylogenies/t__SGB17248_withLongumSubs2/treeshrink/RAxML_bestTree.t__SGB17248.TreeShrink.tre") -> Tree
Tree$tip.label[  grepl("GCA", Tree$tip.label)   ] -> Refs
keep.tip(Tree, c(Refs, longum_infant$sample_id)) -> Tree
longum_infant %>% mutate(Reference = F) %>% rbind(tibble(sample_id = Refs, Country=NA, infant=NA, offset_val=NA,phylo_effect_median=NA, Continent=NA, Reference = T  )) -> longum_infant
longum_infant %>% mutate(infant = ifelse(Reference == T, "Reference", infant) ) -> longum_infant
#There is a repeated ref. I believe it was an issue while fetching in NCBI programaticaly. Need to be remvoed
#longum_infant %>% filter(! sample_name %in% "GCA_001870705") -> longum_infant
#drop.tip(Tree, "GCA_001870705" ) -> Tree

#Get exact infant age
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> Phenotypes
Phenotypes %>% select(ID_anal, infant_age, Age, study_name) %>% rename(sample_id = ID_anal) -> Pheno_infants
left_join(longum_infant, Pheno_infants) -> longum_infant

#get info references
read_tsv("/mnt/project/Samples_tree/SGB/t__SGB17248/2_x/Table_references.tsv") %>% filter(ref %in% Refs) %>% select(ref, subsp) %>% rename(sample_id = ref , Subspecies = subsp) ->  Refs_info
left_join(longum_infant, Refs_info) -> longum_infant
longum_infant %>% mutate( infant = ifelse(Reference == T, NA, infant ) ) -> longum_infant


longum_infant %>% filter(infant == T | Reference == T | Continent == "Africa" ) -> Enrichment_keep 
longum_infant %>% filter(! sample_id %in% Enrichment_keep$sample_id) -> To_sample
sample(To_sample$sample_id,1000) -> Keep2

keep.tip(Tree, c(Keep2, Enrichment_keep$sample_id) ) -> Tree_sub


longest_branch <- which.max(node.depth.edgelength(Tree_sub))
# Root the tree at the midpoint of the longest branch
rooted_tree <- root(Tree_sub, outgroup = Tree_sub$tip.label[longest_branch], resolve.root = T )

#finding the node for infantis
filter(longum_infant,Subspecies == "subsp. suis")$sample_id -> suisref
filter(longum_infant,Subspecies == "subsp. infantis")$sample_id -> infantisref
common_ancestor_infant <- getMRCA(Tree_sub, infantisref[5:10] ) #outgroup infantis
common_ancestor_infant2 <- getMRCA(Tree_sub, suisref ) #outgroup suis
Adults = filter(longum_infant, infant == F)$sample_id
common_ancestor_adult <- getMRCA(Tree_sub, sample(Adults[Adults %in% Tree_sub$tip.label ], 200)  ) 
common_ancestor_adult <- getMRCA(Tree_sub, sample(Adults[Adults %in% Tree_sub$tip.label ], 200)  ) 


rooted_tree$tip.label

ggtree(Tree_sub, layout = "fan",  open.angle=15, size=0.1 ) %<+% longum_infant -> p
p + geom_fruit( geom="geom_tile", mapping = aes(fill=infant), width=0.03,offset=0.1  ) + scale_fill_manual(values = c("TRUE" = "#FF7F7F", "FALSE"="grey"  ), na.value="white" )  + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1) +  scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Reference), width=0.03,offset=0.1) + scale_fill_manual(values = c("TRUE" = "black", "FALSE"="white"  ), na.value="white" )  -> Overall_tree
  #geom_balance(node=6043, fill="blue", color=NA, alpha=0.3)   #-> Overall_tree
#1803 suis
#1824 infantis, not suis
#2024 adult ; 6043 in the whole tree
AdultClade = Tree$tip.label[getDescendants(Tree, node=5831+182)] ; AdultClade[!is.na(AdultClade)] -> AdultClade
longum_infant %>% mutate(Adult_clade = ifelse( sample_id %in% AdultClade, T, F )) -> longum_infant

longum_infant %>% filter(! Reference==T ) %>% rename( Clade = Adult_clade )  %>% mutate(Clade = ifelse(Clade == FALSE, "Infant dominated", "Adult dominated")) %>%
  ggplot(aes(x=Clade, fill=infant )) + geom_bar() + theme_bw() + scale_fill_manual(values = c("TRUE" = "#FF7F7F", "FALSE"="grey"  )) + ylab("Number of samples in clade") + scale_fill_manual(values = c("TRUE" = "#FF7F7F", "FALSE" = "grey"),
  labels = c("TRUE" = "Infant (<=2 years old)", "FALSE" = "Not infant (>2 years old)")) +
  ylab("Number of samples in clade") +  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5, color = "black") + scale_y_log10() -> Proportion_infants_clades
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/InfantCladeBifido_barplot.pdf", Proportion_infants_clades, width = 6.5, height = 4.5)
       
ggtree(Tree , layout = "fan",  open.angle=15, size=0.1 ) %<+% (longum_infant %>% mutate(infant = ifelse(is.na(infant), "Reference", ifelse(infant==T, "Infant (<=2 years old)",  "Not infant (>1 years old)"  ) ) ) ) + geom_balance(node=6014, fill="grey", color=NA, alpha=0.3)-> p
p + geom_tippoint(aes(color=infant), size=0.5)  + scale_color_manual(values = c("#FF7F7F","grey", "black" )) +   new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.1,offset=) + scale_fill_viridis_c(option = "magma") -> TreeBifido
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/TreeBifido_simple.pdf", TreeBifido, width = 6.5, height = 6.5)




longum_infant %>% left_join(. , Metadata %>% select(ID_anal, Age), by=c("sample_id"="ID_anal") ) %>% filter(infant) %>% arrange(desc(Age))

longum_infant %>% filter(infant == T | Reference == T  ) -> Enrichment_keep2 
keep.tip(Tree,  Enrichment_keep2$sample_id ) -> Tree_sub2
ggtree(Tree_sub2 , layout = "fan",  open.angle=15, size=0.1 ) %<+% longum_infant -> p2
p2 + geom_fruit( geom="geom_tile", mapping = aes(fill=infant_age/30), width=0.03,offset=0.1  ) + scale_fill_viridis_c(option = "viridis")  + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1) +  scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Subspecies), width=0.03,offset=0.1) + scale_fill_manual(values = c( "subsp. infantis" = "brown", "subsp. suis"= "darkturquoise", "subsp. longum"="#CAB2D6"  ), na.value="white" ) +
  geom_balance(node= 754+ 735 , fill="red", color=NA, alpha=0.3) + #suis 754 + 3
  geom_balance(node= 754 + 3  , fill="blue", color=NA, alpha=0.3) + #+ #infantes754 + 22 
  geom_balance(node= 754 + 196  , fill="green", color=NA, alpha=0.3) + #754 + 196 
  geom_balance(node= 754 + 182  , fill="grey", color=NA, alpha=0.3) -> infant_tree #754 + 195

drop.tip(Tree, AdultClade ) -> Tree_sub3

#IfnantSpecific = c( getDescendants(Tree_sub2, node = 754 + 3  ), getDescendants(Tree_sub2, node = 754 + 22  ), getDescendants(Tree_sub2, node = 754 + 196 ) )
#Tree_sub2$tip.label[IfnantSpecific] -> IfnantSpecific ; IfnantSpecific[!is.na(IfnantSpecific)] -> IfnantSpecific
#keep.tip(Tree, IfnantSpecific  ) -> Tree_sub3

ggtree(Tree_sub3 , layout = "fan",  open.angle=15, size=0.1 ) %<+% (longum_infant %>% mutate(infant = ifelse( is.na(infant), "Reference", ifelse(infant==T, "Infant (<=2 years old)",   "Not infant (>2 years old)"  ) ) ) )  +
  geom_tippoint(aes(color=infant), size=0.5)  + scale_color_manual(values = c( "#FF7F7F", "grey", "black" )) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=infant_age), width=0.1,offset=0.1  ) + scale_fill_viridis_c(option = "viridis")  + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.1,offset=0.1) +  scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Subspecies), width=0.1,offset=0.1) + scale_fill_manual(values = c( "subsp. infantis" = "brown", "subsp. suis"= "darkturquoise", "subsp. longum"="#CAB2D6"  ), na.value="white" ) +
  geom_balance(node= 399  , fill="#ccfbfe", color=NA, alpha=0.3) + #suis 216+3
  geom_balance(node= 216 + 3  , fill= "#a04668", color=NA, alpha=0.3) -> infant_tree2 #+ #+ #infantes 216 + 22
  #geom_balance(node= 216 + 196  , fill= "#CDD6DD", color=NA, alpha=0.3) -> infant_tree2

#There is an infantis and a longum clustering with the "suis clade"
longum_infant %>% filter(sample_id %in% Tree_sub3$tip.label[ getDescendants(Tree_sub3,399) ] ) %>%
  filter(Reference == T) -> Wrong_assig


ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/BifidoLongum_Infant_all.pdf", Overall_tree)
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/BifidoLongum_Infant_infantonly.pdf", infant_tree)
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/BifidoLongum_Infant_infantonly2.pdf", infant_tree2, width = 8, height=8)  


#enrichment of European samples in adult clade
longum_infant %>% filter(infant==T) %>% group_by(Adult_clade, Continent) %>% summarise(Count = n())  %>% drop_na() %>% ungroup() %>%
  pivot_wider(names_from = Continent, values_from = Count, values_fill = list(Count = 0)) %>%
  column_to_rownames("Adult_clade") %>% as.matrix() %>% fisher.test()


#Check small clade of clonal branches
getDescendants(Tree_sub2,node=754+196) -> Check
#they are not the same samples
Phenotypes %>% filter(ID_anal %in%  Tree_sub2$tip.label[Check] ) %>% select(study_name, sample_id, subject_id, ID_anal) %>% group_by(study_name) %>% summarise(N = n()) %>% arrange(desc(N))
#most samples from Shao, they are not the same baby 
read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/ShaoY.tsv") -> ShaoMeta
read_tsv("/mnt/project/Make_Associations/Phenotypes/Specific/ShaoY_original.tsv") -> ShaoMeta2
ShaoMeta$subject_id %>% sapply(function(x){ str_replace(x, "SID", "") %>% str_replace(., "_ba", "") %>% str_replace(., "_mo", "")   } ) -> ShaoMeta$subject_id
ShaoMeta2$Individual %>%  sapply(function(x){ ifelse(grepl("_", x), str_split(x, "_")[[1]][1], x ) }) -> ShaoMeta2$subject_id
ShaoMeta2 %>% select(subject_id, Time_point, Hospital, Days_in_hospital, Abx_CS_prophylactic, Abx_mother_prior_birth, Abx_mother_labour_IAP, Abx_mother_after_hospital, Abx_Baby_in_hospital, Abx_Baby_after_discharge, Feeding_method) %>% distinct(subject_id, .keep_all=T) %>% left_join(ShaoMeta, .) -> ShaoMeta

ShaoMeta %>% mutate(Rare_clade = ifelse(sample_id %in% Tree_sub2$tip.label[Check] , T, ifelse(sample_id %in% Tree_sub2$tip.label, F, NA) )  ) -> ShaoMeta
ShaoMeta %>% filter(sample_id %in% Tree_sub2$tip.label) %>% filter(! (is.na(Rare_clade) | is.na(infant_age) ) ) %>% glm( as.factor(Rare_clade) ~ born_method + infant_age + Hospital  , . , family=binomial() ) %>% summary()
ShaoMeta %>% filter(sample_id %in% Tree_sub2$tip.label) %>% filter(! (is.na(Rare_clade) | is.na(infant_age) | is.na(Feeding_method) ) ) %>% glm( as.factor(Rare_clade) ~ born_method + infant_age + Feeding_method + Hospital , . , family=binomial() ) %>% summary()

#less likely for vaginal born, slight negative assocation with infant age, 

#enrichment of suis in african
getDescendants(Tree_sub2,node=754+3) -> Check
Phenotypes %>% filter(!is.na(sample_id)) %>% mutate(suis_clade = ifelse(sample_id %in% Tree_sub2$tip.label[Check] , T, ifelse(sample_id %in% Tree_sub2$tip.label, F, NA) )  ) %>% filter(! (is.na(suis_clade) ) ) %>%
  glm( as.factor(suis_clade) ~ Continent , . , family=binomial() ) %>% summary()

#Thre are a couple of genomes (infantis and longum) that fall within the clade of an infantis
#Check if with other methodologies this is also the case

subspecies_colors <- c("subsp. infantis" = "brown", "subsp. suis"= "darkturquoise", "subsp. longum"="#CAB2D6")
clade_colors <- c("suis clade"="#ccfbfe", "infantis clade" = "#a04668", "longum clade" = "grey" )
my_palette <- c("lightgrey", "steelblue")


#Approach A: check ANIs
read_tsv("/mnt/project/Make_Associations/Functional_enrichment/Bifido_Ref/Pangenome/ANI",col_names = F) %>%
  mutate(X1= str_replace(X1, ".gff", ""),X2= str_replace(X2, ".gff", "") ) %>% select(X1, X2, X3) %>%
  spread(X2, X3) -> df_ani
df_ani = df_ani %>% filter(! X1 %in% c("GCA_000155415", "GCA_000196575") ) %>% select( - c("GCA_000155415", "GCA_000196575" ) )
df_ani %>% as.data.frame() %>% column_to_rownames("X1") %>% as.matrix() -> df_ani2
#df_ani2[upper.tri(df_ani2,diag = F)] = df_ani2[lower.tri(df_ani2,diag = F)]
Anno = mutate(df_ani, Clade= ifelse(X1 %in% Wrong_assig$sample_id, 1, 0 )) %>% left_join(longum_infant %>% select(sample_id, Subspecies), by=c("X1" = "sample_id" ) ) %>% select(X1,Clade,Subspecies) %>% as.data.frame() %>% column_to_rownames("X1")
df_ani2[lower.tri(df_ani2)] <- t(df_ani2)[lower.tri(df_ani2)]

#pheatmap::pheatmap(df_ani2, annotation_row=Anno %>% rename(Clade_B.longum_suis = Clade), annotation_colors = list(Subspecies = subspecies_colors, Clade = clade_colors)   ) 
#pheatmap::pheatmap(df_ani2, annotation_row=Anno %>% as.data.frame() %>% column_to_rownames('X1')  , annotation_colors = list(Subspecies = subspecies_colors, Clade = clade_colors)   ) 

#Check the ones of interest
#longum_infant %>% filter(sample_id %in% Tree_sub3$tip.label[ getDescendants(Tree_sub3,216 + 3) ] ) %>%
#  filter(Reference == T) -> Wrong_assig

Anno = mutate(df_ani, Clade= ifelse(X1 %in% Wrong_assig$sample_id, 1, 0 )) %>% left_join(longum_infant %>% select(sample_id, Subspecies), by=c("X1" = "sample_id" ) ) %>% select(X1,Clade,Subspecies)# %>% as.data.frame() %>% column_to_rownames("X1")
Anno %>% mutate(Clade = ifelse(Clade == 1, "suis clade", ifelse(X1 %in%AdultClade, "longum clade", "infantis clade"))) -> Anno
Anno %>% rename(Subs. = Subspecies) -> Anno
pheatmap::pheatmap(df_ani2, annotation_row=Anno %>% as.data.frame() %>% column_to_rownames('X1')  , annotation_colors = list(`Subs.` = subspecies_colors, Clade = clade_colors)   ) 


df_ani %>% filter(X1 %in% Wrong_assig$sample_id ) %>% select(c("X1", one_of(Wrong_assig$sample_id)))
#ANIs > 96%

#Calculate Average ANIs between clades and missannotated taxa
Find_Average_ANI_with_Clade = function(Anno, df_ani2, Genome_ID ){
#  Genome_ID = c("GCA_000092325", "GCA_00155415", "GCA_00196575")

  ANI_df = tibble()
  for (Subclade in unique(Anno$Clade) ){
    Anno %>% filter(Clade == Subclade) -> SC
    df_ani2[Genome_ID, SC$X1] -> ANIs
    Avg_Ani = mean(ANIs)
    ANI_df %>% rbind(tibble(Clade = Subclade, avgANI=Avg_Ani ) ) -> ANI_df
  }
  return(ANI_df)
}
Find_Average_ANI_with_Clade(Anno, df_ani2 ,"GCA_000092325")



#And if taking other infantis references?
longum_infant %>%  filter(Reference == T, grepl("infa", Subspecies ) ) -> infan
c(infan$sample_id, "GCA_00155415", "GCA_00196575" ) -> infant_and_suis
df_ani %>% filter(X1 %in% infant_and_suis ) %>% select(c("X1",one_of(infant_and_suis))) %>% filter(X1 %in%  Wrong_assig$sample_id) %>%
  select(-one_of(Wrong_assig$sample_id))

#Check pangenome

read_csv("/mnt/project/Make_Associations/Functional_enrichment/Bifido_Ref/Pangenome/gene_presence_absence.csv") -> Genes
colnames(Genes)[grepl("GCA", colnames(Genes))] -> Samples 
Genes %>% select(Gene, Samples) -> PangenomeInfo
PangenomeInfo %>% as.data.frame() %>% column_to_rownames("Gene") ->PangenomeInfo2
PangenomeInfo2 %>% mutate_all(~ ifelse(is.na(.), 0, 1)) -> PangenomeInfo2

Anno = mutate(df_ani, Clade= ifelse(X1 %in% Wrong_assig$sample_id, 1, 0 )) %>% left_join(longum_infant %>% select(sample_id, Subspecies), by=c("X1" = "sample_id" ) ) %>% select(X1,Clade,Subspecies)# %>% as.data.frame() %>% column_to_rownames("X1")

Anno %>% rbind( data.frame(X1 = "GCA_900445755", Subspecies="subsp. infantis", Clade=F) ) -> Anno
Anno %>% mutate(Clade = ifelse(Clade == 1, "suis clade", ifelse(X1 %in%AdultClade, "longum clade", "infantis clade"))) -> Anno
Anno %>% rename(Subs. = Subspecies) -> Anno
PangenomeInfo2 %>% select(one_of(Tree_sub3$tip.label, filter(longum_infant, Subspecies == "subsp. longum" )$sample_id )) %>%
  t()  %>% as.data.frame() %>% pheatmap::pheatmap(
    show_colnames = FALSE,annotation_row=Anno %>% as.data.frame() %>% column_to_rownames("X1"), cluster_cols = F,
    annotation_colors = list(Subs. = subspecies_colors, Clade = clade_colors), 
    color = my_palette,legend = FALSE, fontsize_row = 4, 
    annotation_legend = F ) -> Heat
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/Bifido_Pangenome.tiff", Heat, width = 4, height = 3 )
######################################
###Prepare supplementary tables######
######################################

#Sup Table 1: Dataset
##ST 1.1 Participants, accesion, study
##ST 1.2 Summary statistics phenos
##ST 1.3 Trees, number tips, IQTree/RAxML, additonal (with food, with refs, with mice), prevalence SGB per continent, included /not included,. taxonomy
##ST 1.4 Additional genomes:  Food references, Mice datasets, B.longum refs

read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> Phenotypes
Phenotypes %>% select(ID_anal, subject_id, study_name, Country, Continent) -> Phenotypes1.1 #Missing NCBI info
Phenotypes1.1 %>% mutate(study_name = ifelse(study_name=="DAG3", "GacesaR_2022", ifelse(study_name %in% c("LLD", "LLD2") , "AlexandraA_2016", ifelse(study_name=="300OB", "KurilshikovA_2019", ifelse(study_name=="IBD", "VichVilaAV_2018", ifelse(study_name=="300TZFG", "StrazarM_2021", ifelse(study_name=="500FG_FSK", "SchirmerM_2016_FSK",study_name)  ))))) ) -> Phenotypes1.1
#Annonymize dutch IDs
set.seed(812)
To_annonimize = c("GacesaR_2022", "AlexandraA_2016")
Samples_to_annonim = filter(Phenotypes1.1, study_name %in% To_annonimize)$ID_anal
Counter = 1
Annonim_tibble = tibble()
for (A in  sample(Samples_to_annonim, length(Samples_to_annonim)) ){
  An = paste0("Sample_Lifelines_", Counter)
  Counter = Counter + 1
  rbind(Annonim_tibble, tibble(ID_anal=A, Annonym_ID=An) ) -> Annonim_tibble
}
write_tsv(Annonim_tibble, "/mnt/project/Make_Associations/Phenotypes/Dutch_annonym.tsv")
Phenotypes1.1 %>% left_join(Annonim_tibble) %>% mutate(ID_anal = ifelse( !is.na(Annonym_ID), Annonym_ID, ID_anal  ) ) -> Phenotypes1.1
Phenotypes1.1 %>% select(-Annonym_ID) -> Phenotypes1.1


#include NCBI accession
tibble(Accession_number = c("EGAS00001005027", "EGAD00001001991", "EGAS00001003508", "EGAD00001004194") , study_name = c("GacesaR_2022", "AlexandraA_2016","KurilshikovA_2019","VichVilaAV_2018") ) -> Accession_EGA
Phenotypes1.1 %>% mutate( Accession = ifelse(study_name %in% c("GacesaR_2022", "AlexandraA_2016","KurilshikovA_2019","VichVilaAV_2018"), "EGA", "NCBI" )  ) %>% left_join(Accession_EGA) -> Phenotypes1.1
#Get accession for all public cohorts
read_tsv("/mnt/project/Make_Associations/Phenotypes/NCBI_accesions.tsv") %>% rename(ID_anal = id_used, Accession_number=NCBI_accession) -> NCBIs
NCBIs %>% drop_na() %>% left_join(Phenotypes1.1, ., by= c("ID_anal", "study_name"), suffix=c("","_2") ) %>%  mutate(Accession_number = coalesce(Accession_number_2, Accession_number)) %>% select(-Accession_number_2) -> Phenotypes1.1

Phenotypes1.1 %>% mutate(Country = ifelse(study_name == "StrazarM_202", "TZA", Country ) ) %>% mutate(Continent = ifelse(study_name == "StrazarM_202", "Africa", Continent ) ) %>% mutate(study_name = ifelse(study_name == "StrazarM_202", "StrazarM_2021", study_name ) ) -> Phenotypes1.1
Phenotypes1.1 %>% mutate(Accession_number = ifelse(study_name == "StrazarM_2021", "PRJNA686265", Accession_number) ) -> Phenotypes1.1

write_tsv(Phenotypes1.1, "/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table1.1.tsv")

#Summary stats
read_tsv( "/mnt/project/Make_Associations/Association/Results/Info_Cohorts/Prevalence_by_study.tsv") %>%
  mutate(study_name = ifelse(study_name=="DAG3", "GacesaR_2022", ifelse(study_name %in% c("LLD", "LLD2") , "AlexandraA_2016", ifelse(study_name=="300OB", "KurilshikovA_2019", ifelse(study_name=="IBD", "VichVilaAV_2018", ifelse(study_name=="300TZN", "StrazarM_2021", ifelse(study_name=="500FG_FSK", "SchirmerM_2016_FSK",study_name)  ))))) ) -> Individuals_with_pheno
Individuals_with_pheno

read_tsv( "/mnt/project/Make_Associations/Association/Results/Info_Cohorts/SummaryStats_by_study.tsv") %>%
  mutate(study_name = ifelse(study_name=="DAG3", "GacesaR_2022", ifelse(study_name %in% c("LLD", "LLD2") , "AlexandraA_2016", ifelse(study_name=="300OB", "KurilshikovA_2019", ifelse(study_name=="IBD", "VichVilaAV_2018", ifelse(study_name=="300TZN", "StrazarM_2021", ifelse(study_name=="500FG_FSK", "SchirmerM_2016_FSK",study_name)  ))))) ) -> Stats

Individuals_with_pheno$Phenotype %>% unique() -> S2
Stats$Phenotype %>% unique() -> S1

Individuals_with_pheno = Individuals_with_pheno %>% select(Phenotype, Value, study_name, N, Number_individuals_with_pheno, Prevalence) %>% mutate(Phenotype_type = "Categorical")
Stats = Stats %>% select(Phenotype, study_name, Number_individuals_with_pheno,  Mean, SD, Max, Min) %>% mutate(Phenotype_type = "Continuous")

full_join(Individuals_with_pheno, Stats, by=c("Phenotype", "study_name", "Number_individuals_with_pheno", "Phenotype_type") ) -> Pheno_summary

#add category
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypic_categories.tsv") %>% rename(Phenotype = Phenos) %>%
  left_join(Pheno_summary , . ) ->  Pheno_summary

read_tsv("/mnt/project/Make_Associations/Phenotypes/Remove_variables.tsv") -> RM
RM %>% filter(! Remove %in% c("UlcerativeColitis" , "CrohnsDisease" ) ) -> RM
Pheno_summary %>% filter(! Phenotype %in% RM$Remove ) -> Pheno_summary



Pheno_summary %>% write_tsv("/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table1.2.tsv")


#SGBs, N_samples, prevalence, trees and included/not included
read_tsv("/mnt/project/Make_Associations/Association/Tree_sizes2.txt") -> Tree_sizes 
read_tsv("/mnt/project/Make_Associations/Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2", col_names = F) %>% rename(SGB=X1, All_taxonomy=X2) ->  Taxonomy 
Taxonomy %>% separate(All_taxonomy, into = paste0("taxonomy_", 1:7), sep = "\\|") %>% rename(Domain=taxonomy_1, Phylum=taxonomy_2, Class=taxonomy_3, Order=taxonomy_4, Family=taxonomy_5, Genus=taxonomy_6, Species=taxonomy_7  ) -> Taxonomy

str_split(Tree_sizes$Tree, "\\.") %>% lapply(function(x){ x[[2]]} ) %>% unlist() -> SGB_names
str_split(Tree_sizes$Tree, "\\.") %>% lapply(function(x){ x[[1]]} ) %>% unlist() -> Tree_type

Tree_sizes %>% mutate(SGB = str_replace(SGB_names, "t__", ""), Tree_software = str_replace(Tree_type, "_bestTree", "")  )  %>% select(SGB, Tree_software, N_tips) -> Tree_sizes

read_tsv("/mnt/project/SGB_abundances/Result/Prevalence.tsv") -> InfoPrevalence
sapply(InfoPrevalence$SGB, function(x){ str_split(x, "\\|")[[1]] -> y ;  str_replace(y[length(y)], "t__", "")  } ) %>% as.vector() -> SGB_names
InfoPrevalence %>% mutate(SGB = SGB_names) %>% spread(Continent, Prevalence) -> InfoPrevalence
InfoPrevalence[is.na(InfoPrevalence)]=0 
str_split(Taxonomy$Species, ",") %>% lapply(function(x){ if(length(x) > 1 ){ x= x[1] } ; return(x)  } ) %>% unlist() -> Taxonomy$Species
Taxonomy %>% left_join(Tree_sizes) %>% left_join(InfoPrevalence) %>% filter(!is.na(Tree_software)) %>% mutate(Included_analysis = ifelse(N_tips>=300, T, F) ) %>% arrange(desc(N_tips)) -> Taxonomy
Taxonomy %>% rename(prevalence_All=All, prevalence_Africa=Africa, prevalence_Asia=Asia, prevalence_Europe=Europe, prevalence_North_America=North_America, prevalence_Oceania=Oceania, prevalence_South_America=South_America) %>% 
  select(SGB,Domain, Phylum,Class,Order,Family,Genus,Species,Tree_software,N_tips, prevalence_All, prevalence_Africa, prevalence_Asia, prevalence_Europe, prevalence_North_America, prevalence_Oceania, prevalence_South_America, Included_analysis ) -> Taxonomy

Taxonomy %>% write_tsv("/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table1.3.tsv")



#Check which food MAGs were included

#Check which references
read_tsv("/mnt/project/Samples_tree/SGB/t__SGB17248/2_x/Table_references.tsv") -> references_longum
tree_longum = read.tree("/mnt/project/Phylogenies/t__SGB17248_withLongumSubs/treeshrink/RAxML_bestTree.t__SGB17248.TreeShrink.tre")
references_longum %>% filter(ref %in% tree_longum$tip.label ) -> references_tree
references_tree %>% rename(ID_anal = name, NCBI_accession = ref, Taxonomy = subsp  ) %>% mutate( Taxonomy = paste0("B.longum ", Taxonomy ) ) %>% select(-included) -> references_tree
references_tree %>% mutate(study_name= NA, Comments="Isolate assembly") -> references_tree
#Mice datasets
read_tsv("/mnt/project/Make_Associations/Phenotypes/mice_datasets.list", col_names = F) -> mice_data
reference_animalis = read_tsv("/mnt/project/Make_Associations/Phenotypes/Mice/Mice_samples.tsv")
tree_animalis = read.tree("/mnt/project/Phylogenies/t__SGB17278_withFood/treeshrink/RAxML_bestTree.t__SGB17278.TreeShrink.tre")
reference_animalis %>% filter(sample_name %in% tree_animalis$tip.label ) -> references_tree_animalis
references_tree_animalis %>% mutate(NCBI_accession = ifelse(study_name == "BlacherE_2019", "PRJEB32767", ifelse(study_name=="KimMS_2018", NA, "PRJNA540893" ) ) ) -> references_tree_animalis
references_tree_animalis %>% rename(ID_anal = sample_name) %>% mutate( Taxonomy = "SGB17278/Bifidobacterium animalis", Comments="Metagenome" ) -> references_tree_animalis

#Check food
#animalis
#M1231812230     ER_001__ACT__bin.2      ACT     ER      MAG     none    none    SGB17278        GGB10640        FGB3150 60.22   0.0     0.0     http://cmprod1.cibio.unitn.it/databases/MetaRefSGB/resources/Jan21/sequences/00000631/M1231812230.fna.bz2
#ER      none    Environmental   Food production Dairy products  Fermented dairy products        Yoghurt and dietary supplement  none    none    none    none    none    none
#datasetid: ER, MAG
tibble(ID_anal = "M1231812230", study_name="LordanR_2019", Taxonomy="SGB17278/Bifidobacterium animalis",  NCBI_accession = NA, Comments="MetaPhlan4 MAG") -> animalis
#SGB7985 lacoccocus lactis
Process_tree("/mnt/project/Phylogenies/t__SGB7985_withFood/treeshrink/RAxML_bestTree.t__SGB7985.TreeShrink.tre") -> tree_lacto
tree_lacto$tip.label[! tree_lacto$tip.label %in% unique(Phenotypes$ID_anal)] -> references
#M1019106398     BertuzziAS_2018__ERR2212269__bin.4      ERR2212269      BertuzziAS_2018 MAG     none    none    SGB7985 GGB5661 FGB1976 99.25   0.57    0.0     http://cmprod1.cibio.unitn.it/databases/MetaRefSGB/resources/Jan21/sequences/00000631/M1019106398.fna.bz2
tibble(ID_anal = references, study_name="BertuzziAS_2018", NCBI_accession = NA) %>%  mutate(Taxonomy = "SGB7985/Lacoccocus lactis") %>% mutate( Comments="MetaPhlan4 MAG" ) -> lactis


rbind(rbind(rbind(lactis, animalis), references_tree_animalis), references_tree ) -> ST_references

ST_references %>% write_tsv("/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table1.4.tsv")


#Sup Table 2: Batch effect
##ST 2.1: Anpan Schimer results
##ST 2.2: Variability explained PCoA


Anpan_batch = read_csv("/mnt/project/Make_Associations/Association/Results/Summaries/SequencingProtocol_results.csv")
Anpan_batch %>% select(elpd_diff, se_diff, Abs_LowerBound, T_stat, SGB) %>% arrange(elpd_diff) %>% mutate(SGB = str_replace(SGB, "t__", "") ) -> Anpan_batch
write_tsv(Anpan_batch, "/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table2.1.tsv")

read_tsv("/mnt/project/Make_Associations/Association/Results/PCA_partion/Merged_stats.tsv") %>% mutate(SGB = str_replace(SGB, "t__", "") ) %>% select(-File_Name) -> Stats_PCoA
Stats_PCoA %>% write_tsv("/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table2.2.tsv")



#Sup Table 3: Geography
##Geographical effects
##Geographical effects enrichment
##Geographical effects vs traitar / corrected/not corrected
##Geographical effects vs prevalence


Geography_effect = read_tsv("/mnt/project/Make_Associations/Association/Results/Geography_cor/Association_results_merged.tsv") %>% mutate(SGB = str_replace(SGB, "t__", "") ) %>% 
  write_tsv("/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table3.1.tsv")

read_tsv()

Geography_effect_traiter = read_tsv("/mnt/project/Make_Associations/Genome_characteristcs/Results/Traitar_summary_stats.tsv")
write_tsv(Geography_effect_traiter, "/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table3.3.tsv")


Geography_effect_prevalence = read_tsv("/mnt/project/Make_Associations/Genome_characteristcs/Results/EnvironmentPrevalence.tsv")
write_tsv(Geography_effect_prevalence, "/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table3.4.tsv")


#Sup Table 4: Phenotypes
##All summary stats
##Geography statified summary stats (with column indicating country)
##Enrichment
##Geography statified summary stats, with intercept (with column indicating country)
##cohort specific?
read_tsv("/mnt/project/Make_Associations/Phenotypes/Remove_variables.tsv") -> To_remove

AllSummary = read_csv("/mnt/project/Make_Associations/Association/All_results.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  )
AllSummary_f = read_csv("/mnt/project/Make_Associations/Association/All_results_filtered.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  )
AllSummary %>% mutate(Statistical_support = ifelse(Association %in% AllSummary_f$Association, T, F) ) %>% arrange(elpd_diff) %>% filter(!Phenotype %in% To_remove$Remove ) -> AllSummary

write_tsv(AllSummary, "/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table4.1.tsv")



Check_number = function(string){  grepl("^\\d+\\.?\\d*$", string) }


read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Europe_filtered.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) -> sign_Europe
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Asia_filtered.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) -> sign_Asia
#read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica_filtered.csv")%>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  )  -> sign_NorthAmerica
read_csv("/mnt/project/Make_Associations/Association/Results_NA_noCountry_filtered.csv")%>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  )  -> sign_NorthAmerica
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_SouthAmerica_filtered.csv")%>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  )  -> sign_SouthAmerica
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Africa_filtered.csv")%>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  )  -> sign_Africa
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Oceania_filtered.csv")%>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  )  -> sign_Oceania


read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Europe.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) %>% mutate(Statistical_support = ifelse(Association %in% sign_Europe$Association, T, F)) %>% mutate(Continent="Europe") %>% filter(!Check_number(Phenotype)) -> All_Europe
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Asia.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) %>% mutate(Statistical_support = ifelse(Association %in% sign_Asia$Association, T, F))  %>% mutate(Continent="Asia") %>% filter(!Check_number(Phenotype))    -> All_Asia
#read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) %>% mutate(Statistical_support = ifelse(Association %in% sign_NorthAmerica$Association, T, F))  %>% mutate(Continent="North_America") %>% filter(!Check_number(Phenotype))   -> All_NorthAmerica
read_csv("/mnt/project/Make_Associations/Association/Results_NA_noCountry.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) %>% mutate(Statistical_support = ifelse(Association %in% sign_NorthAmerica$Association, T, F))  %>% mutate(Continent="North_America") %>% filter(!Check_number(Phenotype))   -> All_NorthAmerica
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_SouthAmerica.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) %>% mutate(Statistical_support = ifelse(Association %in% sign_SouthAmerica$Association, T, F))  %>% mutate(Continent="SouthAmerica") %>% filter(!Check_number(Phenotype))   -> All_SouthAmerica
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Africa.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) %>% mutate(Statistical_support = ifelse(Association %in% sign_Africa$Association, T, F))  %>% mutate(Continent="Africa") %>% filter(!Check_number(Phenotype))   -> All_Africa
read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Oceania.csv") %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) %>% mutate(Statistical_support = ifelse(Association %in% sign_Oceania$Association, T, F))  %>% mutate(Continent="Oceania")  %>% filter(!Check_number(Phenotype))  -> All_Oceania
                                                                                                                                                                                                                                                                                                                                                                                                                                  
rbind(All_Europe, All_Asia, rbind(All_Africa, rbind(All_NorthAmerica, rbind(All_SouthAmerica)))) %>% arrange(elpd_diff) %>% drop_na() -> StatsperContinent                                                                                                                                                                                                     
write_tsv(StatsperContinent, "/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table4.2.tsv")

                                                                                                                                                                                                               
#Enrichment

Results_enrichment_p= read_tsv("/mnt/project/Make_Associations/Association/Results/Summaries/Enrichment_taxonomiclevel_associations.tsv")
write_tsv(Results_enrichment_p,  "/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table4.3.tsv")




Prepare_associations = function(Path_filtered, Cont = NA){
  Path_unfiltered = str_replace(Path_filtered, "_filtered", "" )
  read_csv(Path_filtered) %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) -> sign
  read_csv(Path_unfiltered) %>% mutate(SGB = str_replace(SGB, "t__", "")  , Association = paste0(SGB, "-", Phenotype ), .before=1  ) %>% mutate(Statistical_support = ifelse(Association %in% sign$Association, T, F)) -> All
  if (! is.na(Cont)) { All %>% mutate(Continent = Cont ) -> All }
  return(arrange(All, elpd_diff) )
}

Prepare_associations("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_Europe_filtered.csv", Cont = "Europe") -> EuropeOffset
Prepare_associations("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_Asia_filtered.csv", Cont = "Asia") -> AsiaOffset
Prepare_associations("/mnt/project/Make_Associations/Association/Results/Summaries/Disease_specific/Results_disease_NorthAmerica_filtered.csv", Cont = "North_America") -> NAOffset
rbind(EuropeOffset, AsiaOffset, rbind(NAOffset) )  %>% arrange(elpd_diff) -> StatsperContinentOffset                                                                                                                                                                                                    
write_tsv(StatsperContinentOffset, "/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table4.4.tsv")




#Sup Table 5: Functional results
##Functional R. ganuvs
##Functional ?
##Functional collisella













read.vcf("/mnt/project/Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Generate_VCF/VCF/t__SGB10068.vcf") -> VCF
library("hierfstat")
library("adegenet")



seq.snp = fasta2DNAbin("/mnt/project/Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Generate_VCF/MSA/t__SGB10068.fa", snpOnly=T, chunkSize = 100)
obj = DNAbin2genind(seq.snp)
meta <- read_tsv('/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv') %>% select(ID_anal, Country)


#Order /filter, #Filter pop. with few samples
meta %>% filter(ID_anal %in% indNames(obj)) -> meta_ob
meta_ob %>% group_by(Country) %>% summarise(N=n()) %>% filter(N > 20) -> To_keep
meta_ob %>% filter(Country %in% To_keep$Country) -> meta_ob
meta_ob %>% distinct(ID_anal, .keep_all = T) -> meta_ob

obj[indNames(obj) %in% meta_ob$ID_anal ] -> obj_f
meta_ob[ match(meta_ob$ID_anal, indNames(obj_f)),  ] -> meta_ob

obj_f$pop = as.factor(meta_ob$Country)

# Convert genind object into hierfstat object
obj.hf = genind2hierfstat(obj_f)

# Various estimation of Fst
genet.dist(obj.hf, method='Nei87', diploid=F) -> FstDistanceMatrix

pairwise.neifst2(obj.hf, diploid=F) -> FstDistanceMatrix
#negative values are consequence of higher intra-population than inter-population differences and should be treated as 0. This can be a result of uneven sample size
FstDistanceMatrix[FstDistanceMatrix<0] = 0


pairwise.neifst2 = function(dat, diploid=F){
  if (is.genind(dat)) { dat <- genind2hierfstat(dat) }
  dat <- dat[order(dat[, 1]), ]
  pops <- unique(dat[, 1])
  Name_pops = pops
  npop <- length(pops)
  fstmat <- matrix(nrow = npop, ncol = npop, dimnames = list(pops,  pops))

  if (is.factor(dat[, 1])) {
    dat[, 1] <- as.numeric(dat[, 1])
    pops <- as.numeric(pops)
  }
  
  for (a in 2:npop) {
    for (b in 1:(a - 1)) {
      print(paste0(Name_pops[a], " vs ",  Name_pops[b]))
      subdat <- dat[dat[, 1] == pops[a] | dat[, 1] == pops[b], ]
      fstmat[a, b] <- fstmat[b, a] <- basic.stats(subdat, diploid = diploid)$overall[8]
    }
  }
  diag(fstmat) <- 0
  return(fstmat)
}    
  

"/mnt/project/Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Fst_matrices/t__SGB10115.csv" %>% read.csv(., header = TRUE) %>% column_to_rownames("X") ->  fst_matrix
intersect(rownames(fst_matrix), rownames(Human_matrix) ) -> countries_check
Human_matrix[countries_check,countries_check] -> Human_matrix
fst_matrix[countries_check,countries_check] -> Bug_matrix
mantel( Human_matrix  ,p,permutations = Per, method="spearman") -> Result_mantel



#Check Fst between continents
file_list <- list.files(path = "/mnt/project/Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Fst_matrices_Nei87/", pattern = "Continent.csv", full.names = TRUE)
Fst_distances = tibble()
for (file_path in file_list) {
  str_split(file_path, "/")[[1]] -> SGBn
  SGBn[[length(SGBn)]] %>% str_replace("_Continent.csv", "") -> SGBn
  # Read the file or perform any desired operation
  data <- read.csv(file_path)
  data %>% as_tibble() %>% rename(Continent1 = X) %>% gather( Continent2, Fst, 2:dim(data)[2] , factor_key=TRUE) %>%
    mutate(SGB = SGBn ) %>% mutate(Continent2 = as.character(Continent2)) -> Fst_info
  IDs = c()
  for (i in seq(1, nrow(Fst_info) ) ){
    paste(sort( c(Fst_info$Continent1[i], Fst_info$Continent2[i]) ), collapse="-") -> i2
    c(IDs, i2) -> IDs
  }
  Fst_info %>% mutate( ID_comparison = IDs) %>% distinct(ID_comparison, .keep_all = T) %>%
    rbind(Fst_distances,. ) -> Fst_distances
  
}


Fst_distances %>% filter(Continent1 != Continent2) %>% filter(!Continent1 %in% c("1", "dumpop") ) %>% mutate(Continent2 = as.character(Continent2)) -> Fst_distances
Fst_distances %>% ggplot(aes(x=ID_comparison, y=Fst)) + geom_boxplot() +  theme_bw() + coord_flip()
Fst_distances %>% filter(!grepl("Oceania", ID_comparison)) %>% filter(!grepl("South_America", ID_comparison)) %>%
  ggplot(aes(x=ID_comparison, y=Fst)) + geom_boxplot() +  theme_bw() + coord_flip()


Fst_distances %>% filter(!grepl("Oceania", ID_comparison)) %>% filter(!grepl("South_America", ID_comparison)) %>% group_by(ID_comparison) %>% summarise(Median = median(Fst), Mean = mean(Fst), SD = sd(Fst) )
Fst_distances %>% filter(!grepl("Oceania", ID_comparison)) %>% filter(!grepl("South_America", ID_comparison)) -> ForTestFst

ForTestFst = list("Africa-Asia"=filter(ForTestFst, ID_comparison=="Africa-Asia")$Fst, "Africa-Europe"=filter(ForTestFst, ID_comparison=="Africa-Europe")$Fst,  "Africa-North_America"=filter(ForTestFst, ID_comparison=="Africa-North_America")$Fst,
     "Asia-Europe"=filter(ForTestFst, ID_comparison=="Asia-Europe")$Fst,  "Asia-North_America"=filter(ForTestFst, ID_comparison=="Asia-North_America")$Fst,  "Europe-North_America"=filter(ForTestFst, ID_comparison=="Europe-North_America")$Fst ) 
dunn.test(ForTestFst)     
    
read_tsv("/mnt/project/Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Summaries_Fst_Mantel_Nei87.tsv") %>% filter(Permutations > 500) %>% select(-BH_FDR) %>% drop_na() %>% mutate(FDR = p.adjust(P, "fdr")) %>% arrange(P)

#Check Fst differences between countries
file_list <- list.files(path = "/mnt/project/Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Fst_matrices_Nei87/", pattern = ".csv", full.names = TRUE)
file_list[!grepl( "Continent" , file_list)] ->file_list_country
Fst_distances_country = tibble()
for (file_path in file_list_country) {
  str_split(file_path, "/")[[1]] -> SGBn
  SGBn[[length(SGBn)]] %>% str_replace(".csv", "") -> SGBn
  # Read the file or perform any desired operation
  data <- read.csv(file_path)
  data %>% as_tibble() %>% rename(Country1 = X) %>% gather( Country2, Fst, 2:dim(data)[2] , factor_key=TRUE) %>%
    mutate(SGB = SGBn ) %>% mutate(Country2 = as.character(Country2)) -> Fst_info
  IDs = c()
  for (i in seq(1, nrow(Fst_info) ) ){
    paste(sort( c(Fst_info$Country1[i], Fst_info$Country2[i]) ), collapse="-") -> i2
    c(IDs, i2) -> IDs
  }
  Fst_info %>% mutate( ID_comparison = IDs) %>% distinct(ID_comparison, .keep_all = T) %>%
    rbind(Fst_distances_country,. ) -> Fst_distances_country
  
}

Fst_distances_country %>% filter(Country1 != Country2) %>% filter(!Country1 %in% c("1", "dumpop") ) %>% mutate(Country2 = as.character(Country2)) -> Fst_distances_country
Fst_distances_country %>% arrange(desc(Fst))
Fst_distances_country %>% group_by(SGB) %>% summarise(M = median(Fst), N = n()) %>% arrange(desc(M)) %>% mutate(SGB = str_replace(SGB, "t__", "") ) -> Median_Fst
read_tsv("/mnt/project/Make_Associations/Association/Results/SupplementaryTables/Table3.1.tsv") -> GeographyEffects

left_join(Median_Fst, GeographyEffects) -> geography_joint
cor.test(geography_joint$M, geography_joint$Rho)

Process_tree("/mnt/project/Symlink_phylo/RAxML_bestTree.t__SGB7865.TreeShrink.tre") -> Tree
Phenotypes %>% filter(Country %in% unique(c(filter(Fst_distances_country,SGB == "t__SGB7865")$Country1, filter(Fst_distances_country,SGB == "t__SGB7865")$Country2 ))) %>% filter(ID_anal %in% Tree$tip.label ) -> Phenos_check
Phenos_check %>% distinct(subject_id, .keep_all = T) -> Phenos_check
keep.tip(Tree, Phenos_check$ID_anal) -> Tree
Phenos_check %>% select(ID_anal,subject_id, Country, Continent, study_name, Age) -> Phenos_check

ggtree(Tree, layout="fan", open.angle=15, size=0.1) %<+% Phenos_check -> p
p +  geom_tiplab(aes(col=Continent, label=Country ), size=1.4) + scale_color_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30") )


read_tsv("/mnt/project/Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Fst_matrices/Human_fst.mat") -> Human_max


#Check if number of polymorphisms is associated with sequencing technology
#for that, we will focus on European sample s(Expecting similar number of polymorphisms) and will associate by year or seq. techn. if available
library(lmerTest)
Clean_names =  function(Name){
  'Samples from 500FG FSK need everything after _ to be removed. If _metaphlan4 is attached at the end of the name, remove'
  if ( grepl("HV", Name)) {
    Name = str_split(Name, "_")[[1]][1]
  }   
  Name = str_replace(Name, "_metaphlan4", "") 
  return(Name)
}
Remove_repeated2 = function(Metadata){
  'Check which samples are repeated (either duplicates or longitudinally) and only picks one from the repeated'
  set.seed(89777)
  Metadata %>% group_by(subject_id) %>% summarise(N = n()) %>% arrange(desc(N)) %>% filter(N > 1 ) -> Replicates
  remove = tibble()
  for (i in Replicates$subject_id){
    Metadata %>% filter(subject_id == i) -> Filtered ; Samples = Filtered$ID_anal
    Keep = sample(Samples, 1)
    remove = c(remove, Samples[!Samples==Keep] ) 
  }
  remove = unlist(remove)
  Metadata %>% filter(! ID_anal %in% remove) -> Metadata
  return(Metadata)
}

read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> Metadata
read_tsv("/mnt/project/Make_Associations/Association/Tree_sizes2.txt") -> Size
Size %>% filter(N_tips >= 300) -> Size
Size$Tree %>% sapply(function(x){str_split(x, "\\.")[[1]][2] }) -> Size$SGB

Location = "/mnt/project/Phylogenies/"
SGBs = list.files(Location, pattern = "*")
Mutation_df = tibble()
for (SGB in SGBs){
  Poly = paste0(Location, SGB, "/", SGB, ".polymorphic") 
  if (! file.exists(Poly)) { next }
  if (! SGB %in% Size$SGB){ next }
  print(SGB)
  read_tsv(Poly) -> Info_poly
  Info_poly$sample %>% sapply(Clean_names) %>% as.vector() -> Info_poly$sample
  
  Info_poly %>% left_join(Metadata, by = c("sample"="ID_anal") ) -> Info_poly
  
  Info_poly %>% filter(Continent == "Europe") -> Info_poly
  Info_poly$study_name %>% sapply(function(x){ str_split(x, "_")[[1]][2] } ) -> Info_poly$Year
  Remove_repeated2(Info_poly %>% rename(ID_anal = sample)) -> Info_poly
  if (length(unique(Info_poly$study_name)) < 5 ){ next }
  if (dim(Info_poly)[1] < 100 ){ next }
  Info_poly %>% mutate(Year = as.numeric(Year)) %>% lmer(percentage_of_polymorphic_sites ~ as.factor(Year) + sequencing_platform + (1|study_name), . ) -> Model2
  
  summary(Model2)$coefficients %>% as.data.frame() %>% rownames_to_column("Feature") %>% mutate(SGB = SGB) %>% rbind(Mutation_df, . ) -> Mutation_df
  
}

as_tibble(Mutation_df) %>% filter(! Feature == "(Intercept)" ) %>% mutate(FDR = p.adjust(`Pr(>|t|)`)) %>% arrange(`Pr(>|t|)`) -> Mutation_df
write_tsv(Mutation_df, "/mnt/project/Make_Associations/Association/Results/Year_and_platform_in_mutationrate.tsv")



#Bristol stool assocation: Eubacterium_rectale
Info_bristol = Process_RDS("/mnt/project/Make_Associations/Association/Results/Continent_stratified/Age,Country_Europe/Models/t__SGB4933_group/BristolType/Model.rds")
Info_bristol %>% filter(BristolType != 0) -> Info_bristol
Tree_Erectale = Process_tree("/mnt/project/Symlink_phylo/IQtree.t__SGB4933_group.TreeShrink.tre")
keep.tip(Tree_Erectale,Info_bristol$sample_id) -> Tree_Erectale
Info_bristol %>% mutate(BristolType_cat = factor(round(BristolType))) -> Info_bristol


GetNodeAncesterClade = function(Tree, df =Info_bristol, Cutoff=0.12, Above= T  ){
  #Get the node that descendents in the individuals with higher phylogenetic effect
  if (Above == T){
    Names = filter(df, phylo_effect_median>Cutoff)$sample_id
  } else {
    Names = filter(df, phylo_effect_median<Cutoff)$sample_id
  }
  match( Names , Tree$tip.label) -> CHECK
  common_ancestor <- getMRCA(Tree, CHECK)
  return(common_ancestor)
}

Info_bristol %>% ggplot(aes(x=BristolType_cat, y=phylo_effect_median )) + geom_boxplot() + theme_bw() + 
  geom_sina(alpha=0.5) 
Info_bristol %>% arrange(desc(phylo_effect_median))

Clade_high = GetNodeAncesterClade(Tree_Erectale, Info_bristol, Cutoff=0.121)
Clade_low = GetNodeAncesterClade(Tree_Erectale, Info_bristol, Cutoff=-0.2)



ggtree(Tree_Erectale, layout="fan", open.angle=15, size=0.1) %<+% Info_bristol -> p
p + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Country), width=0.03,offset=0.1 ) -> p2
p2 + scale_fill_manual(values=c("NLD" = "orange", "ITA"="#1f77b4" )) +
  new_scale_fill()+ geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=BristolType_cat), width=0.03,offset=0.1 ) + 
  scale_fill_manual(values=c("1" = "#ded7c3", "2"="#BEB7A4",  "3"= "#8F8A7B", "4"= "#777367", "5"= "#5F5C52","6"="#302E29", "7"="#000000" )) +
  geom_balance(node=Clade_low, fill="tomato", color=NA, alpha=0.3)


##F. nucleatum
FNuc = Process_tree("/mnt/project/Symlink_phylo//RAxML_bestTree.t__SGB6007.TreeShrink.tre")
Phenos %>% filter(ID_anal %in% FNuc$tip.label ) -> Samples_with_Fuso
Samples_with_Fuso %>% select(study_name, ID_anal, Country, Continent, CRC) %>%
  filter(!is.na(CRC)) -> Samples_with_Fuso





read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_Asia.csv") %>% mutate(Continent = "Asia") %>% rbind( read_csv("/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica.csv") %>% mutate(Continent = "NA") ) -> Compare

#Overall less signal
Compare %>% ggplot(aes(x=Continent, y=log10( abs(elpd_diff)) ) ) + geom_boxplot() + theme_bw()  + geom_hline(yintercept = log10(4) )
#Maybe it is because i am controlling for country...?



#Check the proportion of significant assocaitons between continents
Check_number = function(string){  grepl("^\\d+\\.?\\d*$", string) }


read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_Asia.csv') %>% filter(!Check_number(Phenotype)) -> Asia
read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_Asia_filtered.csv') -> Asia_f


read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica.csv') %>% filter(!Check_number(Phenotype)) -> NorthAm
read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica_filtered.csv') -> NorthAm_f

read_csv('/mnt/project/Make_Associations/Association/Results_NA_noCountry.csv') %>% filter(!Check_number(Phenotype)) -> NA2
read_csv('/mnt/project/Make_Associations/Association/Results_NA_noCountry_filtered.csv') -> NA2_f



read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_Europe.csv') %>% filter(!Check_number(Phenotype)) -> Europe
read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_Europe_filtered.csv') -> Europe_f


dim(Asia_f%>% filter(elpd_diff < -8))[1]  / dim(Asia)[1] 
dim(NorthAm_f%>% filter(elpd_diff < -8))[1] / dim(NorthAm)[1]
dim(NA2_f%>% filter(elpd_diff < -8))[1] / dim(NA2)[1]
dim(Europe_f%>% filter(elpd_diff < -8))[1] / dim(Europe)[1]


length(unique(Europe_f$Phenotype))/length(unique(Europe$Phenotype))
length(unique(Asia_f$Phenotype))/length(unique(Asia$Phenotype))
length(unique(NA2_f$Phenotype))/length(unique(NA2$Phenotype))

Compare_proportion = function(Cut_off_elpd, Cut_off_k ){

  PE = length(unique(filter(Europe_f,(elpd_diff <= Cut_off_elpd & Percentage_badk<=Cut_off_k) )$Phenotype))/length(unique(Europe$Phenotype))
  PA = length(unique(filter(Asia_f, (elpd_diff <= Cut_off_elpd & Percentage_badk<=Cut_off_k) )$Phenotype))/length(unique(Asia$Phenotype))
  PN = length(unique(filter(NA2_f, (elpd_diff <= Cut_off_elpd & Percentage_badk<=Cut_off_k) )$Phenotype))/length(unique(NA2$Phenotype))

  tibble(ELPDdiff_Cutoff=Cut_off_elpd,  Parettok_Cutoff = Cut_off_k, Asia= PA, Europe = PE, North_America = PN ) %>% return()
}

Benchmark_results = tibble()
for (Paretto in  seq(0, 1, by = 0.1) ){
  Compare_proportion(-4, Paretto) %>% rbind(Benchmark_results, . ) -> Benchmark_results
}
Benchmark_results %>% gather(Continent,  Proportion, c(Asia, Europe, North_America) ) %>% filter(ELPDdiff_Cutoff== "-4") %>%
  ggplot(aes(x=Parettok_Cutoff, y=Proportion, col=Continent )) + geom_point() + geom_line() + theme_bw() + ylim(0,1) + scale_color_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30") )



for (ELPD in  seq(-4, -20, by = -1) ){
  Compare_proportion(ELPD, 1) %>% rbind(Benchmark_results, . ) -> Benchmark_results
}
Benchmark_results %>% gather(Continent,  Proportion, c(Asia, Europe, North_America) ) %>% filter(Parettok_Cutoff== 1) %>%
ggplot(aes(x=ELPDdiff_Cutoff, y=Proportion, col=Continent )) + geom_point() + geom_line() + theme_bw() + ylim(0,1) + scale_color_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30") )



Benchmark_results <- NULL
for (Paretto in seq(0, 1, by = 0.1)) {
  for (ELPD in seq(-4, -20, by = -1)) {
    # Apply your function and append to Benchmark_results
    Compare_proportion(ELPD, Paretto) %>%
      rbind(Benchmark_results, .) -> Benchmark_results
  }
}

Benchmark_results %>% gather(Continent,  Proportion, c(Asia, Europe, North_America) ) %>% 
  ggplot(aes(x=ELPDdiff_Cutoff, y=Proportion, col=Continent )) + facet_wrap(~Parettok_Cutoff) +
  geom_point() + geom_line() + theme_bw() + ylim(0,1) + scale_color_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30") )

###DOWNASAMPLING
read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_Asia_filtered.csv')%>% mutate(ID = paste0(Phenotype,"--", SGB) ) -> Asia_f
read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_Asia.csv') %>% filter(!Check_number(Phenotype))%>% mutate(ID = paste0(Phenotype,"--", SGB) ) %>% mutate(Supported = ifelse(ID %in% Asia_f$ID, T, F), Continent ="Asia", N = k_good+k_ok+k_bad+k_verybad ) -> Asia

read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica_filtered.csv')%>% mutate(ID = paste0(Phenotype,"--", SGB) ) -> NorthAm_f
read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_NorthAmerica.csv') %>% filter(!Check_number(Phenotype))%>% mutate(ID = paste0(Phenotype,"--", SGB) )%>% mutate(Supported = ifelse(ID %in% NorthAm_f$ID, T, F), Continent ="NA", N = k_good+k_ok+k_bad+k_verybad )  -> NorthAm

read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_Europe_filtered.csv')%>% mutate(ID = paste0(Phenotype,"--", SGB) ) -> Europe_f
read_csv('/mnt/project/Make_Associations/Association/Results_PerContinent_Europe.csv') %>% filter(!Check_number(Phenotype))%>% mutate(ID = paste0(Phenotype,"--", SGB) )%>% mutate(Supported = ifelse(ID %in% Europe_f$ID, T, F), Continent ="Europe", N = k_good+k_ok+k_bad+k_verybad ) -> Europe

read_csv('/mnt/project/Make_Associations/Association/Results_downsampling_Asia_filtered.csv')%>% mutate(ID = paste0(Phenotype,"--", SGB) ) -> Asia_fd
read_csv('/mnt/project/Make_Associations/Association/Results_downsampling_Asia.csv') %>% filter(!Check_number(Phenotype)) %>% mutate(ID = paste0(Phenotype,"--", SGB) ) %>% mutate(Supported = ifelse(ID %in% Asia_fd$ID, T, F), Continent ="Asia", N = k_good+k_ok+k_bad+k_verybad ) -> Asia_d

read_csv('/mnt/project/Make_Associations/Association/Results_downsampling_NA_filtered.csv')%>% mutate(ID = paste0(Phenotype,"--", SGB) ) -> NorthAm_fd
read_csv('/mnt/project/Make_Associations/Association/Results_downsampling_NA.csv') %>% filter(!Check_number(Phenotype))%>% mutate(ID = paste0(Phenotype,"--", SGB) )%>% mutate(Supported = ifelse(ID %in% NorthAm_fd$ID, T, F), Continent ="NA", N = k_good+k_ok+k_bad+k_verybad ) -> NorthAm_d

read_csv('/mnt/project/Make_Associations/Association/Results_downsampling_Europe_filtered.csv')%>% mutate(ID = paste0(Phenotype,"--", SGB) ) -> Europe_fd
read_csv('/mnt/project/Make_Associations/Association/Results_downsampling_Europe.csv') %>% filter(!Check_number(Phenotype))%>% mutate(ID = paste0(Phenotype,"--", SGB) )%>% mutate(Supported = ifelse(ID %in% Europe_fd$ID, T, F), Continent ="Europe", N = k_good+k_ok+k_bad+k_verybad ) -> Europe_d


rbind(Asia, NorthAm) %>% rbind(Europe) -> Prev
rbind(Asia_d, NorthAm_d) %>% rbind(Europe_d) -> Down
Comparison_down = tibble()
for (Ass in Europe_d$ID ){
  Prev %>% filter(ID == Ass ) %>% select(ID, elpd_diff, phylo_median, Supported, Continent, N) -> Assocaition_prev
  Down %>% filter(ID == Ass )%>% select(ID, elpd_diff, phylo_median, Supported, Continent, N)  -> Assocaition_down
  
  left_join(Assocaition_prev, Assocaition_down, by=c("ID", "Continent" ) , suffix=c('', '_downsample') ) %>% rbind(Comparison_down, . ) -> Comparison_down 
}

Comparison_down %>%  group_by(Supported, Continent ) %>% summarise(N = n()) %>% ungroup() %>% 
  tidyr::spread(Supported, N, fill = 0) %>% mutate(Ratio_Supported= `TRUE` / `FALSE`)%>% mutate(Downsample = F) -> Ratios_nodown
Comparison_down %>%  group_by(Supported_downsample, Continent ) %>% summarise(N = n()) %>% ungroup() %>% 
  tidyr::spread(Supported_downsample, N, fill = 0) %>% mutate(Ratio_Supported= `TRUE` / `FALSE`) %>% mutate(Downsample = T) -> Ratios_down
Comparison_down %>% filter(N_downsample > 100) %>%  group_by(Supported_downsample, Continent ) %>% summarise(N = n()) %>% ungroup() %>% 
  tidyr::spread(Supported_downsample, N, fill = 0) %>% mutate(Ratio_Supported= `TRUE` / `FALSE`) %>% mutate(Downsample = T) -> Ratios_down_cutoff
rbind(Ratios_nodown, Ratios_down) %>% ggplot(aes(x=Continent, y=Ratio_Supported)) + geom_bar(stat="identity") +
  facet_wrap(~Downsample) + theme_bw()


Comparison_down %>%  ggplot(aes(x=Continent, y=log10(N) )) + geom_boxplot(outlier.shape = NA) + theme_bw() + geom_sina(alpha=0.2) + ylab('Number of samples phenotype-phylogeny pair (log10)') -> PanelA
rbind(Ratios_nodown, Ratios_down) %>% mutate(Downsample = ifelse(Downsample == T, 'Downsampled', 'Not downsampled' ) ) %>% ggplot(aes(x=Continent, y=Ratio_Supported)) + geom_bar(stat="identity") +
  facet_wrap(~Downsample) + theme_bw() + ylab('Proportion significant associations') -> PanelB
PanelA | PanelB


###Downsampling percentages EU
Prepare_Ass_table = function( File_stem ){
  Prefix = "/mnt/project/Make_Associations/Association/"
  File_stem = paste0(Prefix, File_stem)
  Ass = paste0(File_stem, "_filtered.csv") %>% read_csv() %>% mutate(ID = paste0(Phenotype,"--", SGB) )
  All =  paste0(File_stem, ".csv") %>% read_csv() %>% mutate(ID = paste0(Phenotype,"--", SGB) ) %>% mutate(Supported = ifelse(ID %in% Ass$ID, T, F), Continent ="Europe", N = k_good+k_ok+k_bad+k_verybad )
  return(All)
}
Eu75 =  Prepare_Ass_table("Results_downsampling_Europe_75") %>% mutate(Continent = "Europe")
Eu50 =  Prepare_Ass_table("Results_downsampling_Europe_50") %>% mutate(Continent = "Europe")
Eu25 =  Prepare_Ass_table("Results_downsampling_Europe_25") %>% mutate(Continent = "Europe")
Eu10 = Prepare_Ass_table("Results_downsampling_Europe_10") %>% mutate(Continent = "Europe")

Eu100 = Prepare_Ass_table("Results_PerContinent_Europe") %>% filter(ID %in% Eu75$ID ) %>% mutate(Continent = "Europe")
NA100 = Prepare_Ass_table("Results_PerContinent_NorthAmerica")%>% filter(ID %in% Eu75$ID ) %>% mutate(Continent = "North_America")
Asia100 = Prepare_Ass_table("Results_PerContinent_Asia")%>% filter(ID %in% Eu75$ID )  %>% mutate(Continent = "Asia")

Plot_down = rbind(Eu75, Eu50) %>% rbind(Eu25) %>% rbind(Eu10) %>% rbind(Eu100) %>% rbind(NA100) %>% rbind(Asia100) %>%
  ggplot(aes(x=N, y= log10(abs(elpd_diff)), col=Continent )) + geom_point() + facet_wrap(~SGB, scales="free", nrow = 1) + geom_hline(yintercept = log10(4) ) + 
  geom_errorbar(aes(ymin = log10(pmax(abs(elpd_diff) - 1.96 * se_diff, 1e-2)), ymax = log10(abs(elpd_diff) + 1.96 * se_diff)),  width = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30") ) + theme_bw()

Eu100 %>% rbind(NA100) %>% rbind(Asia100) %>% select(SGB, elpd_diff, Supported, Continent, N ) %>% arrange(SGB) %>%
  mutate(
    elpd_diff_col = paste0("elpd_diff_", Continent),
    Supported_col = paste0("Supported_", Continent),
    N_col = paste0("N_", Continent)
  ) %>% select(-Continent) %>% 
  pivot_wider(
    names_from = c(elpd_diff_col, Supported_col, N_col),
    values_from = c(elpd_diff, Supported, N)
  ) -> Table_done
colnames(Table_done) = c("SGB", "elpd_diff_Europe", "elpd_diff_NorthAmerica", "elpd_diff_Asia", 
                         "Supported_Europe", "Supported_NorthAmerica", "Supported_Asia", 
                         "N_Europe", "N_NorthAmerica", "N_Asia")
ggsave("/mnt/project/Make_Associations/Git/PhylogeneticAssociation/Analyses_revision/Figures/DownsamplingTrend.pdf",Plot_down, width = 15, height=2)

###IQTREE VS RAXML

#Comparison with IQtrees
Loc = '/mnt/project/Phylogenies_iqtree_redo/Sym/'
Comparison_Stats = tibble()
for (Tree in list.files(Loc) ){
  print(Tree)
  Loc_tree = paste0(Loc, Tree)
  IQ_Tree = Process_tree(Loc_tree)
  SGB = str_split(Tree, "\\.")[[1]][1]
  locraxml = paste0('/mnt/project/Symlink_phylo/RAxML_bestTree.',SGB,'.TreeShrink.tre')
  raxml =  Process_tree(locraxml)
  
  Matrix_iq = cophenetic.phylo(IQ_Tree)
  Matrix_rax =  cophenetic.phylo(raxml)
  Common = intersect( rownames(Matrix_iq), rownames(Matrix_rax)  )
  
  Matrix_iq[Common, Common] -> Matrix_iq
  Matrix_rax[Common, Common] -> Matrix_rax
  
  Stat = vegan::mantel( as.dist(Matrix_iq), as.dist(Matrix_rax), method = "pearson", permutations=999, parallel=5 )
  
  Comparison_Stats = Comparison_Stats %>% rbind( tibble( SGB = SGB,  N_samples = length(Common), Rho_pearson = Stat$statistic, P=Stat$signif, Permutatiopns=Stat$permutations  )  )
}

Comparison_Stats  %>% arrange(desc(Rho_spear)) #%>% gt::gt()

