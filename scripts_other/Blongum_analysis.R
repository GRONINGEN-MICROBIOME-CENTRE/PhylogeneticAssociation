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
Process_tree("/mnt/project/Phylogenies/t__SGB17248_withLongumSubs/treeshrink/RAxML_bestTree.t__SGB17248.TreeShrink.tre") -> Tree
Tree$tip.label[  grepl("GCA", Tree$tip.label)   ] -> Refs
keep.tip(Tree, c(Refs, longum_infant$sample_id)) -> Tree
longum_infant %>% mutate(Reference = F) %>% rbind(tibble(sample_id = Refs, Country=NA, infant=NA, offset_val=NA,phylo_effect_median=NA, Continent=NA, Reference = T  )) -> longum_infant

#Get exact infant age
read_tsv("/mnt/project/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") -> Phenotypes
Phenotypes %>% select(ID_anal, infant_age) %>% rename(sample_id = ID_anal) -> Pheno_infants
left_join(longum_infant, Pheno_infants) -> longum_infant

#get info references
read_tsv("/mnt/project/Samples_tree/SGB/t__SGB17248/2_x/Table_references.tsv") %>% filter(ref %in% Refs) %>% select(ref, subsp) %>% rename(sample_id = ref , Subspecies = subsp) ->  Refs_info
left_join(longum_infant, Refs_info) -> longum_infant


longum_infant %>% filter(infant == T | Reference == T | Continent == "Africa" ) -> Enrichment_keep 
longum_infant %>% filter(! sample_id %in% Enrichment_keep$sample_id) -> To_sample
sample(To_sample$sample_id,1000) -> Keep2

keep.tip(Tree, c(Keep2, Enrichment_keep$sample_id) ) -> Tree_sub


ggtree(Tree_sub, layout = "fan",  open.angle=15, size=0.1 ) %<+% longum_infant -> p
p + geom_fruit( geom="geom_tile", mapping = aes(fill=infant), width=0.03,offset=0.1  ) + scale_fill_manual(values = c("TRUE" = "#FF7F7F", "FALSE"="grey"  ), na.value="white" )  + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1) +  scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Reference), width=0.03,offset=0.1) + scale_fill_manual(values = c("TRUE" = "black", "FALSE"="white"  ), na.value="white" ) -> Overall_tree
  

longum_infant %>% filter(infant == T | Reference == T  ) -> Enrichment_keep2 
keep.tip(Tree,  Enrichment_keep2$sample_id ) -> Tree_sub2
ggtree(Tree_sub2 ,  open.angle=15, size=0.1 ) %<+% longum_infant -> p2
p2 + geom_fruit( geom="geom_tile", mapping = aes(fill=infant_age/30), width=0.03,offset=0.1  ) + scale_fill_viridis_c(option = "viridis")  + 
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.03,offset=0.1) +  scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) +
  new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Subspecies), width=0.03,offset=0.1) + scale_fill_manual(values = c( "subsp. infantis" = "brown", "subsp. suis"= "darkturquoise", "subsp. longum"="#CAB2D6"  ), na.value="white" ) -> infant_tree

ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/BifidoLongum_Infant_all.pdf", Overall_tree)
ggsave("/mnt/project/Make_Associations/Association/Results/General_plots/BifidoLongum_Infant_infantonly.pdf", infant_tree)  
