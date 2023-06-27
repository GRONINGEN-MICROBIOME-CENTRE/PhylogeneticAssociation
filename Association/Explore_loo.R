source("Functions.R")

print("Rscript Explore_loo.R PHENOTYPE_NAME SGB_NAME")
args <- commandArgs(trailingOnly = TRUE)
Pheno = args[1]
SGB = args[2]
if (length(args) == 3){
	Covariates = args[3]
	Covariates = str_split( Covariates, ",")[[1]]
} else{
	Covariates = c()
}

##E.g
#Pheno = "Age"
#SGB = "t__SGB15015"

Files = paste0("Results/Models/",SGB, "/", Pheno)
pattern <- "*_pglmm_fit.RDS"

file_list <- list.files(Files, pattern = pattern, full.names = TRUE)


read_rds(file_list) -> DF
DF$summary() -> DF2


pattern2 = "*_inputs.RData"
file_list = list.files(Files, pattern = pattern2, full.names = TRUE)

load(file_list)


DF2 %>% filter( grepl("phylo", variable))  %>% filter(! grepl("std", variable)) %>% mutate(ID_anal = c("-", as.character(model_input$sample_id)))  -> DF3

Studies = read_tsv("../Phenotypes/Phenotypes_merged2.tsv")

DF3 %>% left_join(Studies %>% select( c("ID_anal", "study_name","Country", "Continent", c(Pheno, Covariates) ))) %>% arrange(desc(abs(median))) -> DF4



DF4 %>% print(n=30)
Out = paste0("Results/Models/",SGB, "/", Pheno, "/DF_posteriors.tsv" )
print(paste0("Saving data in: ", Out))
write_tsv(DF4, Out)


DF4 %>% group_by(Continent) %>% summarise(mean(median), min(median), max(median) ) %>% print()

