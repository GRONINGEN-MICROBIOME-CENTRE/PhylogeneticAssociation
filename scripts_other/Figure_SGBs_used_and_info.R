library(tidyverse)
library(patchwork)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

set.seed(2223)

c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")


Info = read_tsv("Tree_sizes2.txt")
SGB_taxonomy = read_tsv("../Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2", col_names = F )

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
read.tree("../Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk") -> SGB_tree
keep.tip(SGB_tree, Info$SGB_name) -> SGB_tree

#5. Get information about: 
#       Does the tree show sequencing protocol differences?
read_csv("SequencingProtocol_results.csv") -> Info_protocol
Info_protocol %>% select(SGB) %>% mutate( Show_protocol_differences = F ) %>%  mutate(SGB = str_remove(SGB, "t__")) -> Info_protocol
Info %>% left_join(Info_protocol) -> Info
#       Does the tree show geographical differences
read_tsv("Results/Geography_cor/Association_results_merged.tsv") -> Info_geography
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
Root = leaf_nodes[3]
SGB_tree %>% root(., outgroup = Root) -> SGB_tree2

Find_node_ancestor = function(Phylum_name = "Firmicutes", subsample=50, SGB_tree=SGB_tree2){
  set.seed(23422)
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

Find_node_ancestor("Actinobacteria", subsample=500) %>% rbind( Find_node_ancestor("Bacteroidetes") ) %>% 
  rbind( Find_node_ancestor("Firmicutes") ) %>% rbind(Find_node_ancestor("Proteobacteria",subsample=500 )) %>%
  rbind( Find_node_ancestor("Verrucomicrobia",subsample=200) ) %>% rbind(Find_node_ancestor("Lentisphaerae") ) %>% rbind( Find_node_ancestor("Candidatus_Melainabacteria") ) -> Annotation_phyla


ggtree(SGB_tree2, layout="fan", open.angle=15, size=0.1)  %<+% (Info2 %>% mutate( Show_protocol_differences = ifelse(is.na(Show_protocol_differences), "NA", Show_protocol_differences ), Geographical_effect=ifelse(is.na(Geographical_effect), "NA", Geographical_effect ) ) )  -> p
p +   geom_balance(node=filter(Annotation_phyla, Phylum_group=="Bacteroidetes")$node  , fill=color_mapping["Bacteroidetes"], color=NA,  alpha=0.3) +  #Bacteroidetes
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Firmicutes")$node  , fill=color_mapping["Firmicutes"], color=NA,   alpha=0.3) + #Firmicutes
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Proteobacteria")$node  , fill=color_mapping["Proteobacteria"], color=NA,  alpha=0.3) + #Proteobacteria
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Verrucomicrobia")$node  , fill=color_mapping["Verrucomicrobia"], color=NA,   alpha=0.3) + #Verrucomicrobia
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Lentisphaerae")$node  , fill=color_mapping["Lentisphaerae"], color=NA,  alpha=0.3) + #Lentisphaerae
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Candidatus_Melainabacteria")$node  , fill=color_mapping["Candidatus_Melainabacteria"],  color=NA,  alpha=0.3) +  #Candidatus_Melainabacteria
  geom_balance(node=filter(Annotation_phyla, Phylum_group=="Actinobacteria")$node  , fill=color_mapping["Actinobacteria"], alpha=0.3, color=NA ) + #Actinobacteria
  #geom_tippoint(aes_string(col = "Phylum"  ), size= 0.7)  + scale_color_manual(values = color_mapping) + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Included_analysis), width=0.03,offset=0.1 ) + scale_fill_manual( values = c( "FALSE"="grey", "TRUE" = "#E83845") ) + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Log10_number_samples), width=0.03,offset=0.1 ) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Taxonomic_group), width=0.03,offset=0.1 ) + scale_fill_manual( values = c("steelblue", "darkgreen") ) + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Show_protocol_differences), width=0.03,offset=0.1 ) + scale_fill_manual( values = c("FALSE"="grey", "NA"="white" ) ) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Geographical_effect), width=0.03,offset=0.1 ) + scale_fill_manual( values = c("FALSE"="grey", "TRUE"="black", "NA"="white" ) ) -> Plot


ggsave("Results/General_plots/Figure2.pdf", Plot)

  
Info %>% ggplot(aes(x= N_tips)) + geom_histogram() + theme_bw() + scale_x_log10() + geom_vline(xintercept = 300, color= "red",linetype = "dashed", linewidth=2) + xlab("Number of samples per species-tree") -> histogram_all_samples
Info %>% filter(N_tips >= 300) %>% ggplot(aes(x= N_tips)) + geom_histogram() + theme_bw() + xlim(299, 23000) + scale_x_continuous(breaks = seq(300, max(Info$N_tips), by = 1500))  + xlab("Number of samples per species-tree") + theme(axis.text.x = element_text(angle = 45, hjust = 1))  -> histogram_included_samples
histogram_all_samples / histogram_included_samples + plot_layout(widths = c(3, 1)) -> Histogram_plot
ggsave("Results/General_plots/Figure2_histogramSamples.pdf", Histogram_plot)

