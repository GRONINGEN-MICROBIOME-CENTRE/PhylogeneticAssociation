library(tidyverse)
library(ape)
library(vegan)
library(geosphere)
library(maps)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
read_tsv("../Phenotypes/Phenotypes_merged_complete.tsv") %>% filter(! study_name == "ThomasAM_2019_c" ) -> Data
unique(Data$Country) -> C_name
c25[1:length(C_name)] -> C_colors
names(C_colors) = C_name


set.seed(7866)
Run_analysis =  function(Tree_file, Threads, Perm =999 ){
	#Get coordinate info
	#world <- map_data("world") %>% as_tibble() %>% filter(!region == "Antarctica")
	cities = get('world.cities') %>% as_tibble() %>% filter(capital != 0) %>% rename(region=country.etc ) %>% rbind( tibble(name="Suva", region="Fiji", pop=NA,lat=-18.1405, long=178.44, capital=1 ))
	read_csv("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Make_Associations/Phenotypes/Countries.csv", col_names = F) %>% rename(Country=X1, region=X2 ) -> Translation
	Translation %>% rbind(tibble(Country = "KOR", region="Korea South" )) -> Translation
	#Get tree
	Clean_names =  function(Name){
  		'Samples from 500FG FSK need everything after _ to be removed. If _metaphlan4 is attached at the end of the name, remove'
  		if ( grepl("HV", Name)) { Name = str_split(Name, "_")[[1]][1] }   
  		Name = str_replace(Name, "_metaphlan4", "") 
  		return(Name)
	}
	read.tree(Tree_file) -> Tree
	if (class(Tree) ==  "multiPhylo"){ Tree = Tree[[1]] }
	Tree = midpoint.root(Tree)
	Tree$tip.label = Tree$tip.label %>% sapply(Clean_names)
	#Get data
	read_tsv("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Make_Associations/Phenotypes/Phenotypes_merged_complete.tsv") %>% filter(! study_name == "ThomasAM_2019_c" ) -> Data
	Data %>% filter(ID_anal %in% Tree$tip.label) -> Data
	#remove repeated samples
	Data %>% distinct(subject_id, .keep_all = T) %>% distinct(ID_anal, .keep_all = T)  -> Data2
	#Merge data and geospatial information
	Data2 %>% select(ID_anal, Country, Continent, study_name) %>% left_join( cities %>% left_join(Translation)  ) %>% group_by(ID_anal, region) %>% sample_n(1) -> Data3

	#Merge data and tree information
	Data3 %>% filter(ID_anal %in% Tree$tip.label) -> Data_anal
	keep.tip(Tree, Data_anal$ID_anal) -> Tree
	order_indices <- match(Tree$tip.label, as.character(Data_anal$ID_anal))
	Data_anal[order_indices,] -> Data_anal

	#Make plot
	ggtree(Tree, layout = "circular")  %<+% Data_anal + geom_tippoint(aes(col=Country), size=0.5 ) + scale_color_manual( values = C_colors  ) + 
    new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes(fill=Continent), width=0.01,offset=0.1 ) + scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30"))
  

	#calculate spatial and phylogenetic distances
	Data_anal %>% ungroup() %>% select(long, lat) %>% as.matrix()  %>% distm(. ,  fun = distVincentySphere) -> distance_matrix
	cophenetic.phylo(Tree) -> Distance_tree
	#perform mantel test
	mantel( as.dist(Distance_tree), as.dist(distance_matrix), method = "pearson", permutations=Perm, parallel=Threads ) -> Test #0.3 r
	tibble(Rho= Test$statistic, P=Test$signif, Permutations= Test$permutations) -> Results

	return(list(Results, Country_plot))

}


SGB_from_tree = function(filename){
	split_filename <- unlist(strsplit(filename, "\\."))
	# Extract the second element
	extracted_string <- split_filename[2]
}

args <- commandArgs(trailingOnly = TRUE)
Results = tibble()
if ( length(args) > 0 ) {
	Tree_file = args[1]
	SGB = SGB_from_tree(Tree_file)
	print(SGB)
	Out = paste0("Results/Geography_cor/", SGB, ".tsv")
	if(file.exists(Out)){ print("Already done") ; q() }
	Run_analysis(Tree_file, Threads =10, Perm=2000 ) -> list_results
	list_results[[1]] %>% mutate(SGB=SGB) -> Results
	write_tsv(Results, Out)
	list_results[[2]] %>% ggsave( str_replace(Out, ".tsv", ".pdf"), . )
} else{
	tre_files <- list.files(path = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Symlink_phylo/", pattern = "\\.tre$", full.names = TRUE)
	for (Tree_file in Files){
		SGB = SGB_from_tree(Tree_file)
		print(SGB)
		Out = paste0("Results/Geography_cor/", SGB, ".tsv")
		if(file.exists(Out)){ print("Already done") ; next }
		Run_analysis(Tree_file, Threads =10, Perm=2000) -> list_results
		list_results[[1]] %>% mutate(SGB=SGB) -> Result
		write_tsv(Result, Out)
		Result %>% rbind(Results, . ) -> Results
		list_results[[2]] %>% ggsave( str_replace(Out, ".tsv", ".pdf"), . )
	}
	write_tsv(Results, "Results/Geography_cor/All_cor.tsv")

}




