library(tidyverse)
library(ape)
library(vegan)
library(phytools)

#####################
#Bifido animalis#####
#####################

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
Distances_to_plot %>%
  ggplot(aes(x=Distance, y =distance, fill=Relationship )) + geom_boxplot() + theme_bw() + scale_y_continuous(trans='log10') + ylab("Phylogenetic distance") + xlab("Distance between") + scale_fill_manual( values= c("Different continent" = "#E31A1C", "Same continent; different country"= "#FDBF6F", "Same country; different study"= "steelblue4", "Same study"= "black"  )  ) + coord_flip() +
  theme(
    text = element_text(size = 16),  # Set the text size for all elements
    plot.title = element_text(size = 20, face = "bold"),  # Set title text size
    axis.title = element_text(size = 18, face = "bold"),  # Set axis title text size
    axis.text = element_text(size = 14),  # Set axis text size
    axis.ticks = element_line(size = 1.5),  # Set axis tick size
    legend.text = element_text(size = 14),  # Set legend text size
    legend.title = element_text(size = 16, face = "bold")  # Set legend title text size
  )

Distances_to_plot %>% mutate( Distance = ifelse(Distance == "Human_to_Human", Relationship, Distance ) ) %>% aov(distance ~ Distance, data = .) -> anova_results
TukeyHSD(anova_results)-> tukey_results
tukey_results$Distance %>% as.data.frame() %>% rownames_to_column("Comparison") %>% as_tibble() %>% arrange(`p adj`) %>% mutate(Sign = ifelse(`p adj`<0.05, T, F ) )
#Everything is sign different than Mice
#Nothing is sign diff than food
#Between countries, only same_study vs Different_continent is not sign. 
#these stats are inflated anyhow. Would need permutations if we need to support with stats
