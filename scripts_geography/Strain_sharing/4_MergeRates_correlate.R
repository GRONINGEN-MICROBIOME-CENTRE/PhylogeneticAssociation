library(tidyverse)


file_list <- list.files(path = "./Transmission/", pattern = "Tranmission_proportions.tsv", recursive = TRUE, full.names = TRUE)

if (! file.exists("Tranmission_rates.tsv")){ 
df_transmission = tibble()
for (file in file_list){
	df = read_csv(file)
	#df %>% select(-`...1`) -> df
	SGB = basename(dirname(file))
	df %>% mutate(SGB = SGB) -> df
	df_transmission %>% rbind(df) -> df_transmission
}
	write_tsv(df_transmission, "Tranmission_rates.tsv")
} else{
	read_tsv("Tranmission_rates.tsv") -> df_transmission
}



Phylo = read_tsv("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Make_Associations/Association/Results/Geography_cor/Association_results_merged.tsv")

#df_transmission %>% mutate(Relationship = ifelse(motherbaby==T, "MotherBaby", ifelse(Family==T, "Family", ifelse(Intrastudy==T, "Intrastudy", "Interstudy")))) -> df_transmission

#Intercountry,Intracountry_diffStudy,Intracountry_sameStudy,Family,motherbaby,Shared
df_transmission %>% mutate(Relationship = ifelse(motherbaby==T, "MotherBaby", ifelse(Family==T, "Family", ifelse(Intercountry==T,"Intercountry", ifelse(  Intracountry_diffStudy == T, "Intracountry_diffStudy", ifelse(Intracountry_sameStudy == T, "Intracountry_sameStudy" , NA ) ) ) ) ) ) -> df_transmission


df_transmission %>% select(Relationship, Shared, SGB) %>% spread(Relationship, Shared) -> df_transmission

print("===Missing SGBS=====")
Phylo %>% filter(!SGB %in% df_transmission$SGB) %>% print()
left_join(df_transmission, Phylo) -> df_transmission


print("===Check the range of shareness ratios===")
print( paste0("Same country: ", max(df_transmission$Intracountry_diffStudy)," - " ,min(df_transmission$Intracountry_diffStudy) , " - ", median(df_transmission$Intracountry_diffStudy) ) )
print( paste0("Diff country: ", max(df_transmission$Intercountry)," - " ,min(df_transmission$Intercountry)," - ",  median(df_transmission$Intercountry) ) )



print("====Testing sharedness vs geographical effect====")
cor.test(df_transmission$Rho, df_transmission$Intracountry_diffStudy) %>% print()
for (i in c("Family","Intracountry_diffStudy","Intracountry_sameStudy","Intercountry","MotherBaby")){
	Model = paste0("Rho ~ ", i)
	df_transmission %>% lm(Model, . ) %>% summary() -> Result
	print(Result$adj.r.squared)
	Result$coefficients %>% as.data.frame() %>% rownames_to_column("Relationship") %>% as_tibble() %>% filter(Relationship == i) %>% print()

}

df_transmission %>% filter(FDR<0.05) -> df_transmission2
cor.test(df_transmission2$Rho, df_transmission2$Intracountry_diffStudy) %>% print()
print("====Testing sharedness vs geographical effect: Repeat only on signifcant SGBs====")
for (i in c("Family","Intracountry_diffStudy","Intracountry_sameStudy","Intercountry","MotherBaby")){
        Model = paste0("Rho ~ ", i)
        df_transmission %>% filter(FDR<0.05)  %>% lm(Model, . ) %>% summary() -> Result
	print(Result$adj.r.squared)
        Result$coefficients %>% as.data.frame() %>% rownames_to_column("Relationship") %>% as_tibble() %>% filter(Relationship == i) %>% print()

}





#Is there an enrichment of sharedness in specific taxonomic level?
SGB_taxonomy = read_tsv("../Make_Associations/Phenotypes/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2", col_names = F ) %>% rename(SGB=X1, Taxonomy=X2)
df_transmission %>% mutate(SGB = str_replace(SGB, "t__", "") ) -> df_transmission
SGB_taxonomy <- SGB_taxonomy %>% separate(Taxonomy, into = c("domain","phylum", "class", "order", "family", "genera", "species"), sep = "\\|", remove = FALSE)
left_join(df_transmission, SGB_taxonomy) ->  df_transmission

library(fgsea)
taxa_list = list()
for (i in c("phylum", "class", "order", "family", "genera") ){
  taxa_list_pre <- split(df_transmission$SGB, df_transmission[i] %>% as_vector() %>% as.vector() )
  taxa_list <- c(taxa_list, taxa_list_pre)
}
df_transmission %>% arrange(desc(Intracountry_diffStudy))-> R1
#Prepare ranks
ranks = R1$Intracountry_diffStudy
names(ranks) = as.character(R1$SGB)
#Run GSEA
fgseaRes <- fgsea::fgsea(taxa_list, ranks)
as_tibble(fgseaRes) %>% arrange(padj) -> Assoc_rho_taxa
print(Assoc_rho_taxa)






#Merge with previous  estimates


#Intercountry	Intracountry_diffStudy	Intracountry_sameStudy	Family	motherbaby	Shared	SGB

Trans =read_tsv("../Make_Associations/Phenotypes/VallesColomerTransmission.tsv")

df_transmission %>% mutate(SGB = str_replace(SGB, "t__", "")) %>%  left_join(Trans) -> All_info


print("====Correlations between sharedness estimates computed and published in Valles-Colomer====")
print("Household")
cor.test( All_info$Family, All_info$SGB_household_transmissibility, na.action="omit" ) %>% print()
print("Intrastudy")
cor.test( All_info$Intracountry_sameStudy, All_info$SGB_intradataset_transmissibility, na.action="omit" ) %>% print()
paste(median(All_info$Intracountry_sameStudy), median(drop_na(All_info)$Intracountry_sameStudy),  median( drop_na(All_info)$SGB_intradataset_transmissibility) )   %>% print()


print("Mother-baby")
cor.test( All_info$MotherBaby, All_info$`SGB_mother-infant_transmissibility`, na.action="omit" ) %>% print()



print("Test with Mireia's")
All_info %>% drop_na() -> All_info_m
print("Household")
cor.test(All_info_m$Rho, All_info_m$SGB_household_transmissibility) %>% print()
print("Intrastudy")
cor.test(All_info_m$Rho, All_info_m$SGB_intradataset_transmissibility) %>% print()
print("Intestudy")
cor.test(All_info_m$Rho, All_info_m$`SGB_mother-infant_transmissibility`) %>% print()
q()


print("===Check differences between Quartiles of interstudy transmission=====")
df_transmission %>% mutate(Quartile_sharedness = ntile(Interstudy, 4)) %>% filter(Quartile_sharedness %in% c("1", "4")) -> df_transmission
df_transmission %>% ggplot(aes(x = as.factor(Quartile_sharedness), y = Rho)) + geom_boxplot() + labs(x = "Quartile of SGB sharedness", y = "Rho Geography vs Tree") + theme_minimal()  -> Plot
df_transmission %>% lm(Rho ~ Quartile_sharedness, . ) %>% summary()

ggsave("Quartiles_sharedness_VS_rho.pdf", Plot)





