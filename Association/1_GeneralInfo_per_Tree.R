library(dplyr)
library(ape)
library(stringr)
library(readr)

Info_tree = tibble()

Location ="../Phylogenies/Symlink_phylo/"
file_list <- list.files(Location, full.names = TRUE)

for (file in file_list) {
	if (file == "../Phylogenies/Symlink_phylo//Make_sym.sh" ){ next }
	print(file)
	#RAxML_bestTree.t__SGB13950.TreeShrink.tre
	Tree = read.tree(file)
	SGB = str_split(file, '\\.')[[1]][4]
	if (class(Tree) == "multiPhylo" ){ Tree = Tree[[1]] } 
	I = tibble(SGB = SGB, Number_leafs = length(Tree$tip.label))
	Info_tree %>% rbind( I ) -> Info_tree
}

Info_tree %>% write_tsv("Results/Tree_stats.tsv")

