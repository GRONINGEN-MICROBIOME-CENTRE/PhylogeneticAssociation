library(tidyverse)



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

To_plot =  c("Alkaline phosphatase", "Gram positive", "Trehalose", "Motile", "L-Rhamnose", "Catalase", "Casein hydrolysis", "Salicin", "Glucose fermenter")
Traitar_merged %>% select( c(Rho, FDR, SGB, To_plot) ) %>% gather(Trait, Available, c(To_plot) ) %>% filter(! Available == 0.5 ) %>% left_join(SGB_taxonomy %>% select(SGB, family) ) %>% mutate(Lachnospiraceae = ifelse(family == "f__Lachnospiraceae", T, F) ) %>%
  mutate(Available = ifelse(Available == 0, "No", "Yes")  ) %>%
  ggplot(aes(y=Rho, x= as.factor(Available), col=Lachnospiraceae  )) + geom_boxplot(outlier.shape = NA) + geom_sina(alpha=0.5) + theme_bw() + facet_wrap(~Trait ) + scale_color_manual(values= c("FALSE" = "blue", "TRUE" = "red" ) ) + xlab("Predicted trait") +
  theme(
    text = element_text(size = 16),  # Set the text size for all elements
    plot.title = element_text(size = 20, face = "bold"),  # Set title text size
    axis.title = element_text(size = 18, face = "bold"),  # Set axis title text size
    axis.text = element_text(size = 14),  # Set axis text size
    axis.ticks = element_line(size = 1.5),  # Set axis tick size
    legend.text = element_text(size = 14),  # Set legend text size
    legend.title = element_text(size = 16, face = "bold")  # Set legend title text size
  )


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
geom_text_repel( aes(label = ifelse(FDR < 0.05, as.character(Annotation), "")), box.padding = 0.5, point.padding = 0.5,  size = 5,   force = 5,  box.color = NA, direction = "both"  ) 

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


Prev_env %>% filter(environment %in% c("Stool_Ancient") ) %>% ggplot(aes(x= as.factor(prevalence_bin) , y = Rho )) + geom_boxplot(outlier.shape = NA) + geom_sina() + theme_bw() + ylab("Geographical effect (rho)") + xlab("SGB present in Ancient stool samples") +theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))



