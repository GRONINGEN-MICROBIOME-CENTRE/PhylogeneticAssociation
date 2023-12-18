library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  File_name <- args[1]
  # Your code that uses the input_value
  print(paste("Received input:", File_name))
} else {
  File_name = "All_results"
  #File_name = "Geography_results"
}



read_csv(paste0(File_name, ".csv") ) -> DF
#DF %>% mutate(P = 1 - pnorm(T_stat),  FDR = p.adjust(P, "fdr")) -> DF
DF %>% mutate( Percentage_badk = 100* ( (k_bad+k_verybad)/(k_good+k_ok+k_bad+k_verybad))  ) -> DF
#DF %>% filter(SGB %in% c("t__SGB9203","t__SGB14232", "t__SGB14025", "t__SGB15067" ) ) %>% select(SGB, Phenotype, Percentage_badk, Abs_LowerBound , elpd_diff)  %>% print()
#q()

####If want to try with a frequentist approach... Not recommended
Frequentist_filtering = function(x){
	#1. Check significant assocaitons
	DF %>% filter(FDR < 0.05 & elpd_diff < -4 )  %>% arrange(elpd_diff) %>% group_by(Phenotype) %>% summarise(N = n()) %>% arrange( desc(N) ) %>% print(n=36)
	DF %>% filter(FDR < 0.05 & elpd_diff < -4 )  %>% arrange(elpd_diff) %>% group_by(SGB) %>% summarise(N = n()) %>% arrange( desc(N) ) %>% print(n=36)
	#2. Make volcano plot
	DF %>% mutate(Significance = ifelse(FDR>0.05, "Not_significant", ifelse(elpd_diff < -4, "Significant", "Significant_smallVariability" ) ) ) %>% mutate(BadK = ifelse( Percentage_badk ==0 , "0%",  ifelse(Percentage_badk < 1, ">0%, <1%", ">=1%"  )) ) %>%
  	ggplot(aes(x=elpd_diff, y=-log10(P), col= Significance, shape=BadK  ) ) + geom_point() + theme_bw() + geom_vline(xintercept = -4) -> volcano_plot
	ggsave( paste0("Results/General_plots/Assocations_volcano_",File_name,".pdf"), volcano_plot)
	#3. How many bad Ks are there: Seem to be driven by outliers in the phenotype.
	Check_bad = function(DF, threshold){
	  DF %>% filter(FDR < 0.05 & elpd_diff < -4 & Percentage_badk >= threshold  ) -> BadK_assocations
	  BadK_assocations %>% group_by(Phenotype) %>% summarise(N = n()) %>% arrange(desc(N)) %>% print()
	  BadK_assocations %>% group_by(SGB) %>% summarise(N = n()) %>% arrange(desc(N)) %>% print()
  	print(BadK_assocations)
	}
	Check_bad(DF, 1)
	cat("We will use threshold of:\n1. FDR < 0.05, based on the statistic obtained from elpd_diff/se_diff, which with enough samples is believed to follow a standard normal distribution\n2. A elpd_diff of less than -4. This filters out associations where the amount of variability that is explained is rather low. Yes, it is better than the NULL, but the LOO is bad. Might be interesting to play a bit with this with continuous variables... If a phenotype varaibility is influenced by many factors, the observed explained variability might be rather small. Maybe we will explore some of them.\n3. We will not check assocations were there are over 1% of the samples being highly influencial (pareto k bad or very bad). This usually happens when there are outliers in the values.
")

	print("Complete: 1. All assocations, 2.Number Phenotypes, 3. Number Bugs ")
	DF %>% dim()
	DF %>% group_by(Phenotype) %>% summarise(n()) %>% dim()
	DF %>% group_by(SGB) %>% summarise(n()) %>% dim()

	print("FDR: 1. All assocations, 2.Number Phenotypes, 3. Number Bugs")
	DF %>% filter(FDR < 0.05) %>% dim()
	DF %>% filter(FDR < 0.05) %>%  group_by(Phenotype) %>% summarise(n()) %>% dim()
	DF %>% filter(FDR < 0.05) %>% group_by(SGB) %>% summarise(n()) %>% dim()

	print("FDR and elpd_diff: 1. All assocations, 2.Number Phenotypes, 3. Number Bugs")
	DF %>% filter(FDR < 0.05 & elpd_diff < -4) %>% dim()
	DF %>% filter(FDR < 0.05 & elpd_diff < -4) %>%  group_by(Phenotype) %>% summarise(n()) %>% dim()
	DF %>% filter(FDR < 0.05 & elpd_diff < -4) %>% group_by(SGB) %>% summarise(n()) %>% dim()

	print("FDR and elpd_diff and pareto: 1. All assocations, 2.Number Phenotypes, 3. Number Bugs")
	DF %>% filter(FDR < 0.05 & elpd_diff < -4 & Percentage_badk < 1) %>% dim()
	DF %>% filter(FDR < 0.05 & elpd_diff < -4 & Percentage_badk < 1) %>%  group_by(Phenotype) %>% summarise(n()) %>% dim()
	DF %>% filter(FDR < 0.05 & elpd_diff < -4 & Percentage_badk < 1) %>% group_by(SGB) %>% summarise(n()) %>% dim()

	DF %>% filter(FDR < 0.05 & elpd_diff < -4 & Percentage_badk < 1) %>% select(Phenotype, SGB, elpd_diff, Percentage_badk) -> Filtered_df

}


DF %>% mutate( CI95_overlaps0 = ifelse( Abs_LowerBound > 0, F, T) ) %>% ggplot(aes(x=phylo_median, y=elpd_diff, col=CI95_overlaps0, shape=Percentage_badk<0)) + geom_point() + theme_bw() + geom_hline(yintercept = -4) -> associations_summary_plot
ggsave( paste0("Results/General_plots/Assocations_summary_",File_name,".pdf"), associations_summary_plot)

print("Complete: 1. All assocations, 2.Number Phenotypes, 3. Number Bugs ")
DF %>% dim()
DF %>% group_by(Phenotype) %>% summarise(n()) %>% dim()
DF %>% group_by(SGB) %>% summarise(n()) %>% dim()

print("Abs_LowerBound > 0 : 1. All assocations, 2.Number Phenotypes, 3. Number Bugs")
DF %>% filter(Abs_LowerBound > 0) %>% dim()
DF %>% filter(Abs_LowerBound > 0) %>%  group_by(Phenotype) %>% summarise(n()) %>% dim()
DF %>% filter(Abs_LowerBound > 0) %>% group_by(SGB) %>% summarise(n()) %>% dim()

print("Abs_LowerBound > 0 and elpd_diff: 1. All assocations, 2.Number Phenotypes, 3. Number Bugs")
DF %>% filter(Abs_LowerBound > 0 & elpd_diff < -4) %>% dim()
DF %>% filter(Abs_LowerBound > 0 & elpd_diff < -4) %>%  group_by(Phenotype) %>% summarise(n()) %>% dim()
DF %>% filter(Abs_LowerBound > 0 & elpd_diff < -4) %>% group_by(SGB) %>% summarise(n()) %>% dim()


print("Abs_LowerBound > 0 and elpd_diff and pareto: 1. All assocations, 2.Number Phenotypes, 3. Number Bugs")
DF %>% filter(Abs_LowerBound > 0 & elpd_diff < -4 & Percentage_badk < 1) %>% dim()
DF %>% filter(Abs_LowerBound > 0 & elpd_diff < -4 & Percentage_badk < 1) %>%  group_by(Phenotype) %>% summarise(n()) %>% dim()
DF %>% filter(Abs_LowerBound > 0 & elpd_diff < -4 & Percentage_badk < 1) %>% group_by(SGB) %>% summarise(n()) %>% dim()




DF %>% filter(Abs_LowerBound > 0 & elpd_diff < -4 & Percentage_badk < 1) %>% select(Phenotype, SGB, elpd_diff, Percentage_badk,phylo_mean, phylo_sd ) -> Filtered_df
Filtered_df %>% group_by(Phenotype) %>% summarise(N_associations = n() ) %>%  ggplot(aes(x= reorder(Phenotype, N_associations ), y=N_associations )) + geom_bar(stat = "identity") + theme_bw() + coord_flip() + theme(axis.text=element_text(size=10) ) -> Bar_plot
ggsave( paste0("Results/General_plots/",File_name ,"N_assocations.pdf") , Bar_plot)

Filtered_df %>% group_by(Phenotype) %>% summarise(N = n()) %>% arrange(desc(N)) %>% print(n=100)
#Filtered_df  %>% print(n=100)


write_csv(Filtered_df, paste0(File_name, "_filtered.csv") )




if (File_name == "Geography_results"){
Remove_names=str_split("T_stat,k_good,k_ok,k_bad,k_verybad,phylo_median,phylo_mad,phylo_q5,phylo_q95,phylo_rhat,phylo_ess_bulk,phylo_ess_tail", ",")[[1]]
wide <- DF %>% select(-c(Remove_names)) %>%
                 pivot_wider(names_from = Phenotype,
                values_from = c(elpd_diff,se_diff,Abs_LowerBound,Percentage_badk,phylo_mean,phylo_sd),
                names_sep = "_")

wide_filterd <- Filtered_df %>% 
                 pivot_wider(names_from = Phenotype,
                values_from = c(elpd_diff,Percentage_badk,phylo_mean,phylo_sd),
                names_sep = "_")

print("Wide filtered")
print(wide_filterd)
print(wide_filterd %>% dim())
print("who is missing?")
wide %>% filter(! SGB %in% wide_filterd$SGB)  %>% write_tsv("Results/Geography_plots/Not_selected.tsv")


#Based on European
print("Wide complete")
wide %>% dim() %>% print()
print("Wide removal boundary takes 0")
wide %>% filter(Abs_LowerBound_Europe > 0) %>% dim() %>% print()
print("Wide removal boundary and elpd_diff")
wide %>% filter(Abs_LowerBound_Europe > 0 & elpd_diff_Europe < -4) %>% dim() %>% print()
print("Wide removal boundary, elpd_diff and bad paretoK")
wide %>% filter(Abs_LowerBound_Europe > 0 & elpd_diff_Europe < -4 & Percentage_badk_Europe < 1) %>% dim() %>% print()


}










