library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(ape)

#anpan
library(cmdstanr)
set_cmdstan_path("/groups/umcg-gastrocol/tmp01/Sergio/.cmdstan/cmdstan-2.32.2/")
library(anpan)
#phylr
library(ggtree)
library(phyr)


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
  if (Fam == "binomial"){ Meta[Phenotype_name] =  ifelse(as.factor(as_vector(Metadata[Phenotype_name])) == "1", T, F) } #If it is not logical, it gives me error
  
  Labels = Meta$ID_anal[!is.na(Meta[Phenotype_name] ) ]
  Meta %>% filter(ID_anal %in% Labels) -> Meta
  keep.tip(Tree, Meta$ID_anal) -> Tree

  Meta %>% select( c("ID_anal", Phenotype_name, COV) ) %>% rename( sample_id=ID_anal ) -> Meta

  print("Fitting PGLMM")
  anpan_pglmm(Meta, Tree, outcome = Phenotype_name, covariates = COV, omit_na = T, family=Fam, out_dir=Model_save, save_object =T,parallel_chains = 4 ) -> result
  write_rds(result, paste0(Model_save, "Model.rds"))
  #print(result$loo$comparison)

  print("Saving plot")
  plot = plot_tree_with_post(tree_file  = Tree, meta_file  = Meta, fit   = result$pglmm_fit, outcome  = Phenotype_name, labels=result$model_input$sample_id, covariates=COV, omit_na =T )
  write_rds(plot, str_replace(Plot, ".pdf", ".rds"))
  ggsave(Plot, plot)


  as_tibble(result$loo$comparison)[2,] %>% mutate(Phenotype = Phenotype_name, covariates= paste(COV, collapse="+")  ) -> Result
  return(Result)

}


Function_fit_phyr = function(Phenotype_name, Tree, Meta, Plot , Fam="gaussian", COV = c(1), ID_random=NULL ){
  if (Fam == "binomial"){ Meta[Phenotype_name] =  ifelse(as.factor(as_vector(Metadata[Phenotype_name])) == "1", T, F) } #If it is not logical, it gives me error

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

