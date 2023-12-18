library(tidyverse)
library(hierfstat)
library(adegenet)
library(vegan)
library(ape)

cat("Conda env: anpan_env ")
cat("Usage: FASTA_to_FST.R SGB_ID")


Location = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/"
args <- commandArgs(trailingOnly = TRUE)
SGB = args[1]
if (length(args) == 2){
	Grouping = args[2]
} else {
Grouping = "Country"
}


FASTA = paste0(Location,  "Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Generate_VCF/MSA/", SGB, ".fa" )
Phenotype = paste0(Location,"Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv")
#Outputs
if (Grouping == "Continent"){
Bug_fst_file = paste0(Location, "Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Fst_matrices_Nei87/",SGB,"_Continent.csv")
} else {
Bug_fst_file = paste0(Location, "Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Fst_matrices_Nei87/",SGB,".csv")
}
Location_human = paste0(Location, "Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Fst_matrices/Human_fst.mat")
Association_summary = paste0(Location, "Make_Associations/HumanPhylogeny_vs_BacteriaPhylogeny/Association_Fst_Nei87/", SGB, ".tsv")

print(paste0("Input: ", FASTA) )
print(paste0("Output: ", Bug_fst_file) )


#Read data

seq.snp = fasta2DNAbin(FASTA, snpOnly=T, chunkSize = 100)
obj = DNAbin2genind(seq.snp)
meta <- read_tsv(Phenotype) %>% select( c("ID_anal", Grouping) )


#Order /filter, #Filter pop. with few samples
meta %>% filter(ID_anal %in% indNames(obj)) -> meta_ob
meta_ob %>% group_by(!!sym(Grouping)) %>% summarise(N=n()) %>% filter(N > 20) -> To_keep
print("Samples to use: "); print(To_keep)

meta_ob %>% filter( !!sym(Grouping) %in% as_vector(To_keep[Grouping]) ) -> meta_ob
meta_ob %>% distinct(ID_anal, .keep_all = T) -> meta_ob

obj[indNames(obj) %in% meta_ob$ID_anal ] -> obj_f
meta_ob[ match(meta_ob$ID_anal, indNames(obj_f)),  ] -> meta_ob

obj_f$pop = as.factor(as_vector(meta_ob[Grouping]))


# Convert genind object into hierfstat object
obj.hf = genind2hierfstat(obj_f)


#usually we would do this method for computing pairwise Fst, howhere the distance matrix is odd, so we will directly use a slightly changed neifst function
#genet.dist(obj.hf, method='Nei87', diploid=F) -> FstDistanceMatrix

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
      print(paste0(Name_pops[a], " vs ", Name_pops[b]))
      subdat <- dat[dat[, 1] == pops[a] | dat[, 1] == pops[b], ]
      fstmat[a, b] <- fstmat[b, a] <- basic.stats(subdat, diploid = diploid)$overall[8]
    }
  }
  diag(fstmat) <- 0
  #negative values are consequence of higher intra-population than inter-population differences and should be treated as 0. This can be a result of uneven sample size
  fstmat[fstmat<0] = 0
  return(fstmat)
}    


if (!file.exists(Bug_fst_file)){
        pairwise.neifst2(obj.hf, diploid=F)  -> fst_matrix
	write.csv(fst_matrix, file = Bug_fst_file, row.names = T)
} else {
        fst_matrix <- read.csv(Bug_fst_file, header = TRUE) %>% column_to_rownames("X")
}

if (Grouping == "Continent"){ q() }

#Compare with human Fst
Human_matrix = read_tsv(Location_human) %>% as.data.frame() %>% column_to_rownames("...1")
intersect(rownames(fst_matrix), rownames(Human_matrix) ) -> countries_check
Human_matrix[countries_check,countries_check] -> Human_matrix
fst_matrix[countries_check,countries_check] -> Bug_matrix
Per = 100000
mantel( Human_matrix  , Bug_matrix,permutations = Per, method="spearman") -> Result_mantel
print(Result_mantel)
print(warnings())
tibble(Rho=Result_mantel$statistic,  P = Result_mantel$signif, Permutations = Result_mantel$permutations, Size_human=length(colnames(Human_matrix)), Size_bug=length(colnames(Bug_matrix))  ) %>% write_tsv(Association_summary)


