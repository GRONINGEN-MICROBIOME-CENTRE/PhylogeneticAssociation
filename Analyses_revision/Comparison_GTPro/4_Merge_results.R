library(tidyverse)

Merged_fi = 'Correlations_per_taxa.tsv'

Merge = function(){
    Files = list.files('Results_mantel', pattern = 'tsv$', full.names=T)
    do.call(rbind, lapply(Files, read_tsv )  ) -> Merged_files
    Merged_files %>% mutate(FDR = p.adjust(P, "fdr") ) -> Merged_files
    Merged_files %>% arrange(desc(Rho)) -> Merged_files
    write_tsv(Merged_files, Merged_fi)
}

Repeat = F
if (file.exists(Merged_fi) & Repeat == F  ){
    Merged_files = read_tsv(Merged_fi)
} else {
    Merged_files = Merge()
}

print('Number of species tested')
dim(Merged_files)[1] %>% print()

##Filter out files with very few N
Merged_files %>%  filter(N >= 50) %>% mutate(FDR = p.adjust(P, "fdr") ) -> Merged_files
print('Number of species tested (over N threshold)')
dim(Merged_files)[1] %>% print()

print('N FDR significant files')
Merged_files %>% group_by(FDR<0.05) %>% summarise(n()) %>% print()
print('Check distribution of number of samples')
Merged_files$N %>% summary() %>% print()
print('Check distribution of rho')
Merged_files$Rho %>% summary() %>% print()
