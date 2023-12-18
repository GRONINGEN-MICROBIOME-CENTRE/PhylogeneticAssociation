library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(ape)

#anpan
library(cmdstanr)
#set_cmdstan_path("/groups/umcg-gastrocol/tmp01/Sergio/.cmdstan/cmdstan-2.32.2/")
library(anpan)
#phylr
#library(ggtree)
#library(phyr)


Compare_shrinked =  function(Tree1, Tree2){
  'Compares two trees, show which branches are in Tree1 that are not in Tree2'	
  print(length(Tree1$tip.label))
  print(length(Tree2$tip.label) )
  if ( length(Tree1$tip.label) !=  length(Tree2$tip.label) ){
    print( paste0("Tree has been shrinked: ",length(Tree1$tip.label), " to ", length(Tree2$tip.label) ) ) 
    print("Missing samples")
    tibble(Tip = Tree1$tip.label) %>% filter(! Tip %in% Tree2$tip.label) %>% print()    
    
  }   
}

Clean_names =  function(Name){
  'Samples from 500FG FSK need everything after _ to be removed. If _metaphlan4 is attached at the end of the name, remove'
  if ( grepl("HV", Name)) {
    Name = str_split(Name, "_")[[1]][1]
  }   
  Name = str_replace(Name, "_metaphlan4", "") 
  return(Name)
}


Remove_repeated = function(Metadata, pruned_tree){
  'Check which samples are repeated (either duplicates or longitudinally) and only picks one from the repeated'
  set.seed(89777)
  Metadata %>% group_by(subject_id) %>% summarise(N = n()) %>% arrange(desc(N)) %>% filter(N > 1 ) -> Replicates
  remove = tibble()
  for (i in Replicates$subject_id){
    Metadata %>% filter(subject_id == i) -> Filtered ; Samples = Filtered$ID_anal
    Keep = sample(Samples, 1)
    remove = c(remove, Samples[!Samples==Keep] ) 
  }
  remove = unlist(remove)
  Metadata %>% filter(! ID_anal %in% remove) -> Metadata
  pruned_tree %>% drop.tip(remove) -> pruned_tree
  return(list(Metadata, pruned_tree))
}



Function_fit_anpan = function(Phenotype_name, Tree, Meta, Plot, Model_save , Fam="gaussian", COV = NULL){
  'Prepares for fit using anpan and saves data'
  #Prepare binomial in a format anpan does not give error
  if (Fam == "binomial"){
	if (class(as_vector(Meta[Phenotype_name])) != "logical"){
	Meta[Phenotype_name] =  ifelse(as.factor(as_vector(Metadata[Phenotype_name])) == "1", T, F) } #If it is not logical, it gives me error
  }
  if (Phenotype_name %in% COV ) { COV = COV[!COV==Phenotype_name]  }
  Labels = Meta$ID_anal[!is.na(Meta[Phenotype_name] ) ]
  Meta %>% filter(ID_anal %in% Labels) -> Meta
  keep.tip(Tree, Meta$ID_anal) -> Tree

  Meta %>% select( c("ID_anal", Phenotype_name, COV) ) %>% rename( sample_id=ID_anal ) -> Meta
  print("Fitting PGLMM")
  
  COV2 = c()
  for (Covariate in COV){
	Meta[Covariate] %>% unique() %>% dim() -> Dim
	if (Dim[1] > 1 ){ COV2 = c(COV2, Covariate) }
  }
  #print(Meta)
  #print( summary( as_vector(Meta[Phenotype] )))
  #q()
  anpan_pglmm(Meta, Tree, outcome = Phenotype_name, covariates = COV2, omit_na = T, family=Fam, out_dir=Model_save, save_object =T,parallel_chains = 4 ) -> result
  write_rds(result, paste0(Model_save, "Model.rds"))
  #print(result$loo$comparison)

  print("Saving plot")
  plot = plot_tree_with_post(tree_file  = Tree, meta_file  = Meta, fit   = result$pglmm_fit, outcome  = Phenotype_name, labels=result$model_input$sample_id, covariates=COV2, omit_na =T )
  write_rds(plot, str_replace(Plot, ".pdf", ".rds"))
  ggsave(Plot, plot)

  print("Saving pareto")
  #Save pareto
  Do_pareto(result, paste0(Model_save, "pareto.tsv"))
  print("Saving phylogenetic estimate")
  Get_phylo_estimate(result, paste0(Model_save, "phylo_estimates.tsv"))

  as_tibble(result$loo$comparison)[2,] %>% mutate(Phenotype = Phenotype_name, covariates= paste(COV2, collapse="+")  ) -> Result
  return(Result)

}


Function_fit_phyr = function(Phenotype_name, Tree, Meta, Plot , Fam="gaussian", COV = c(1), ID_random=NULL ){
  if (Fam == "binomial"){ Meta[Phenotype_name] =  ifelse(as.factor(as_vector(Metadata[Phenotype_name])) == "1", T, F) 
  } else {  Meta[Phenotype_name] = as.numeric(as_vector(Metadata[Phenotype_name]))  } #If it is not logical, it gives me error

  Labels = Meta$ID_anal[!is.na(Meta[Phenotype_name] ) ]
  Meta %>% filter(ID_anal %in% Labels) -> Meta
  keep.tip(Tree, Meta$ID_anal) -> Tree
  
  Meta %>% select( c("ID_anal","subject_id" , Phenotype_name, COV) ) %>% rename( sample_id=ID_anal ) -> Meta

  print("Fitting PGLMM")
  Formula_b = as.formula( paste0(Phenotype_name, " ~ ", paste(c(COV, ID_random), collapse="+")  ))
  Formula_c = as.formula( paste0(Phenotype_name, " ~ ", paste(c(COV, ID_random, "(1|sample_id__)"), collapse="+") ))
  if (is.null(ID_random) == F) {
  	#Base model
	phyr::pglmm(Formula_b , data = Meta, family = Fam, REML = FALSE ) -> BM
  	#Complete model
  	phyr::pglmm(Formula_c , data = Meta, family = Fam, REML = FALSE, cov_ranef = list(sample_id = Tree), s2.init = c(BM$ss)^2 ) -> CM
	BM_lik = BM$logLik
  } else {
	glm(Formula_b, data = Meta, family = Fam )  -> BM 
        phyr::pglmm(Formula_c , data = Meta, family = Fam, REML = FALSE, cov_ranef = list(sample_id = Tree) ) -> CM
	BM_lik = as.numeric(logLik(BM))
  }
  #Comput LRT #https://daijiang.github.io/phyr/reference/pglmm.html : Examples 
  # Twice the log of the likelihood ratio (LLR) between models is asymptotically (as n gets large) chi-square distributed
  pvalue <- pchisq(2*(CM$logLik - BM_lik), df = 1, lower.tail = FALSE)
  tibble(Phenotype = Phenotype_name, P=pvalue, Covariates= paste( c(COV, ID_random), collapse="+") ) %>% mutate(Phenotype = Phenotype_name) -> Result

  print("Saving plot")
  ggtree(Tree,layout = "circular"  ) -> p
  p %<+% Meta + theme(legend.position="bottom") + geom_tippoint(aes_string(col = Phenotype_name  ), size= 0.7) -> p
  ggsave(Plot, p)


  return(Result)

}

Get_phylo_estimate =  function(RDS_model, Output){
        RDS_model$pglmm_fit -> fit
        fit$summary("sigma_phylo") -> Params
        if( is.null(Params) ){ return() }
        write_tsv(Params, Output)
}


Get_pareto = function(vector){
        breaks <- c(-Inf, 0.5, 0.7, 1, Inf)
        labels <- c("good", "ok", "bad", "very bad")
        # Create the table
        result <- table(cut(vector, breaks, labels = labels, include.lowest = TRUE))
        # Calculate the percentage
        percentage <- prop.table(result) * 100
        # Calculate the minimum
        minimum <- tapply(vector, cut(vector, breaks, labels = labels, include.lowest = TRUE), min, na.rm = TRUE)
        # Combine the results into a data frame
        table_df <- data.frame(Count = as.integer(result),
                       Pct. = paste0(format(percentage, digits = 1), "%"),
                       Min. = minimum )
        # Print the table
        return(table_df) }

Get_pareto_wide = function(vector){
        # Define the breaks for categorization
        breaks <- c(-Inf, 0.5, 0.7, 1, Inf)
        # Define the labels for each category
        labels <- c("good", "ok", "bad", "very bad")
        # Create the table
        result <- table(cut(vector, breaks, labels = labels, include.lowest = TRUE))
        # Extract the "Count" column in wide format
        count_wide <- as.integer(result)
        # Print the count in wide format
        return(count_wide)
}

Do_pareto = function(RDS_model, Output){
        #For each rds model
        Loo = RDS_model$loo$pglmm_loo
        Pareto = Loo$diagnostics$pareto_k

        Get_pareto_wide(Pareto) -> pareto_values #"good", "ok", "bad", "very bad" 
        names(pareto_values) = c("good", "ok", "bad", "very bad")
        t(pareto_values) %>% as.data.frame() -> pareto_values
        write_tsv(pareto_values, Output)
	print("Pareto saved")
}


Tree_plot =  function(SGB, Metadata_location="../Phenotypes/Phenotypes_merged_complete.tsv", Pheno="Continent", Show_plot =F, Overwrite=F ){
	library(ape)
	library(tidyverse)
	library(ggtree)
	library(phytools)
	Plot_name = paste0("Results/Geography_plots/", SGB, ".pdf")	
	if (Overwrite == F){
		if ( file.exists(Plot_name) ){ return() }
	}

	print("Reading tree")
	Tree = paste0("../../Symlink_phylo/RAxML_bestTree.", SGB, ".TreeShrink.tre")
	if (! file.exists(Tree) ){ Tree = paste0("../../Symlink_phylo/IQtree.", SGB, ".TreeShrink.tre")
	}
	Tree = read.tree(Tree)
	if (class(Tree) == "multiPhylo" ){ Tree = Tree[[1]] }
	Tree = midpoint.root(Tree)
	#Rename branches
	Tree$tip.label %>% sapply(. , Clean_names) -> New_names
	Tree$tip.label = New_names
	print("Reading meta")
	Metadata = read_tsv(Metadata_location) %>% filter(! study_name == "ThomasAM_2019_c" )
	Metadata %>% group_by(ID_anal) %>% sample_n(1) -> Metadata
	#Trimming metadata
	Metadata %>% filter(ID_anal %in% Tree$tip.label) -> Metadata
	Metadata = Metadata[! is.na(Metadata[Pheno] ), ]


	# Pruning tree
	print("Pruning tree")
	pruned_tree <- keep.tip(Tree, Metadata$ID_anal )

	# Sort Metadata
	print("Matching data")
	Metadata %>% filter(ID_anal %in% pruned_tree$tip.label) -> Metadata
	Metadata %>% group_by(ID_anal) %>% summarise(N = n()) %>% arrange(desc(N))
	Metadata <- arrange(Metadata, match(ID_anal, pruned_tree$tip.label ))
	
	#Removal longitudinal samples
	print("Removing longitudinal samples")
	Unrepeated = Remove_repeated(Metadata, pruned_tree)
   	Metadata = Unrepeated[[1]] ; pruned_tree = Unrepeated[[2]]

    	#Make plot
   	print("Plotting")
	#It is important to have the tip labels in the first column so that they match the tree
    	Metadata %>% select(c("ID_anal", Pheno)) -> Metadata2
	#ggtree plot: Make circular, attache the metadata with the operator %<+% , color by Pheno. Color pallete is given manually.
  	ggtree(pruned_tree,layout = "circular"  ) -> p
  	p %<+% Metadata2 + theme(legend.position="bottom") + geom_tippoint(aes_string(col = Pheno  ), size= 0.7) -> p
  	p + scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#8c564b"
 ),aesthetics = "colour", na.value="white") -> p
  	ggsave(Plot_name, p)
	if (Show_plot == T){ plot(p) }

}

Geography_results_to_wide = function(Filtered = F){
	if (Filtered == F){
	wide <- DF %>% select(-c(se_diff, T_stat)) %>%
 		 pivot_wider(names_from = Phenotype,
              	values_from = c(elpd_diff, Abs_LowerBound, k_good, k_ok, k_bad, k_verybad),
              	names_sep = "_")
	} else {
	 wide <- DF  %>%
                 pivot_wider(names_from = Phenotype,
                values_from = c(elpd_diff, Percentage_badk),
                names_sep = "_")
	}
	return(wide)
}


root.by.outgroup <- function(tree.unrooted) {
  longest.edge <- which.max(tree.unrooted$edge.length)
  long.nodes <- tree.unrooted$edge[longest.edge,]     #this should usually include one tip
  new.outgroup <- long.nodes[long.nodes < Ntip(tree.unrooted)]
  tree.rooted <- root(tree.unrooted, outgroup=new.outgroup, resolve.root=T)
  tree.rooted
}

Tree_plot_v2 =  function(SGB, Metadata_location="../Phenotypes/Phenotypes_merged_complete.tsv", Pheno="Continent", Cov=c("Age", "Continent"), Special_cov=c(),Plot_name=NA,Model_name=NA, Tree=NA  ){                                                          
        library(ape)
        library(tidyverse)
        library(ggtree)
        library(phytools)
        library(ggtreeExtra) #Allows to add metadata to circular plots
        library(ggnewscale) #Allos for multiple fill scales within one plot
        library(viridis) #color palette

	print("---Plotting Tree---")

        print("Reading model")
	if (is.na(Plot_name)){
        	Model_name = paste0("Results/",paste(Cov, collapse=","),"/Models/",SGB,"/",Pheno,"/Model.rds")
       		Model = read_rds(Model_name)
	} else { Model = Model_name }
        print("Getting information from model")
        Model$model_input -> Info
        Model$pglmm_fit$summary() -> Fit
        print("Preparing information table")
        Fit %>% filter( grepl("phylo", variable))  %>% filter(! grepl("std", variable)) %>% mutate(sample_id = factor(c("-", as.character(Info$sample_id)))) %>% select(sample_id, median) %>% rename(phylo_effect_median = median) -> Phylo
        left_join(Info, Phylo) %>% mutate(sample_id = as.character(sample_id)) -> Info
	if (! class(Tree) == "phylo" ) {
 	       print("Reding tree")
		Tree_file = paste0("../../Symlink_phylo/RAxML_bestTree.", SGB, ".TreeShrink.tre")
		if (!file.exists(Tree_file)){ Tree_file = paste0("../../Symlink_phylo/IQtree.",SGB, ".TreeShrink.tre")  }
	       	Tree = read.tree(Tree_file)
		if (class(Tree) == "multiPhylo" ){ Tree = Tree[[1]] }
       		print("Processing tree")
        	Tree = midpoint.root(Tree) #root in midpoint
        	#Tree =root.by.outgroup(Tree) #root in longest branch
        	Tree$tip.label %>% sapply(. , Clean_names) -> New_names
        	Tree$tip.label = New_names
        	Tree= keep.tip(Tree, Info$sample_id)
	}
	#Sometimes not all covariates are included, for intance if there is only one continent. We need to account for this
	Cov2 = c()
  	for (Covariate in Cov){
		if (Covariate %in% colnames(Info) ){ 
			Cov2 = c(Cov2, Covariate)
			if (Covariate == "Age" & Pheno=="Age" ){ Cov2 = c(Cov2, "Age_group") }
		}
	}
	if ("Age" %in% colnames(Info)){
		#Info %>% mutate(Age_group = ifelse(Age ==0, "0", ifelse(Age < 3, "(0-3)", ifelse(Age < 10, "[3-10)", ifelse(Age<20, "[10-20)", ifelse(Age<40, "[20-40)", ifelse(Age<60, "[40-60)", ifelse(Age<80, "[60-80)", ifelse(Age>=80, ">80", Age)))))))  ) ) %>% mutate(Age_group = factor(Age_group, levels=c( "0", "(0-3)", "[3-10)", "[10-20)", "[20-40)", "[40-60)", "[60-80)", ">80" )))   -> Info
		Info %>% mutate(Age_group = ifelse(Age ==0, "<1", ifelse(Age < 15, "(0-15)", ifelse( Age < 30, "[15, 30)", ifelse( Age < 70, "[30, 70)", ">=70") )))) %>% mutate(Age_group = factor(Age_group, levels=c( "<1", "(0-15)", "[15, 30)", "[30, 70)", ">=70"  ))) -> Info
	}
	

        #Plot
	Info = as_tibble(Info)
	Info %>% drop_na() -> Info
        print("Making tree plot")
        ggtree(Tree, layout="fan", open.angle=15, size=0.1) %<+% Info -> p
        p + geom_fruit( geom="geom_tile", mapping = aes(fill=phylo_effect_median), width=0.03,offset=0.1) + scale_fill_viridis_c(option = "magma") -> p2
        for (Ph in unique(c(Cov2, Pheno))){
		#if (Ph == "Age"){ Ph = "Age_group" }

                p2 + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes_string(fill=Ph), width=0.03,offset=0.1 ) -> p2
                if (Ph == "Continent"){
                        p2 + scale_fill_manual(values=c("Europe" = "#1f77b4", "North_America" = "#ff7f0e", "South_America" = "#2ca02c", "Asia" = "#d62728", "Oceania" = "#9467bd", "Africa" = "#FFCE30")) -> p2
                } else if (Ph == "Age"){
                        p2 + scale_fill_viridis_c(option = "viridis") -> p2
		 } else if (Ph == "Age_group"){
			#p2 + scale_fill_manual(values = c("0" = "#e41a1c",  "(0-3)"= "#377eb8", "[3-10)" = "#4daf4a", "[10-20)"= "#984ea3", "[20-40)"="#ff7f00", "[40-60)"="#a65628", "[60-80)"="#f781bf",  ">80"="#999999"), na.value="white" ) -> p2
			p2 + scale_fill_manual(values = c("<1" = "#5E1DB5",  "(0-15)"= "#5A3C82", "[15, 30)" = "#0E27E8", "[30, 70)"= "#EBB649",  ">=70"="#B5731D" ), na.value="white" ) -> p2
                } else if (class( as_vector(Info[Ph]) )  == "logical"){
                        p2 + scale_fill_manual( values = c("grey", "#E83845") ) -> p2
                } else if (class( as_vector(Info[Ph]) )  != "logical"){
                        p2 + scale_colour_gradient2() -> p2
                }
        }

	if (length(Special_cov) > 0){
                read_tsv("../Phenotypes/Phenotypes_merged_complete.tsv") %>% select(-sample_id) %>% rename(sample_id = ID_anal) %>% select( c("sample_id", Special_cov) ) -> To_merge
                left_join(Info, To_merge) -> Info
        
                for (Ph2 in unique(Special_Cov) ){
                        p2 + new_scale_fill() + geom_fruit( geom="geom_tile", mapping = aes_string(fill=Ph), width=0.03,offset=0.1 ) + scale_fill_manual( values = c("grey", "black"), na.value = "white" ) -> p2
                }
        }
        p2 + theme( legend.background=element_rect(fill=NA),legend.title=element_text(size=6), legend.text=element_text(size=4.5), legend.spacing.y = unit(0.02, "cm") ) -> p2
	print("Saving")
	if (is.na(Plot_name)){
        	Plot_name = paste0("Results/General_plots/Posterior_and_covariates/", paste(Cov, collapse=","),"/", paste0(SGB,"_",Pheno) )     
       }
	print(paste0("Saving: ", Plot_name) )
        ggsave( paste0(Plot_name, ".pdf") , p2)
        ggsave( paste0(Plot_name, ".png") , p2)
}



