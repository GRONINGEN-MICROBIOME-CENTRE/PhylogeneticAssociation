library(tidyverse)
library(ape)
library(ggforce)
library(coin)
library(vegan)
library(phytools)

setwd("/scratch/p300317/GlobalHumanGutMicrobiomeGeneticDiversity_2024/Technical_batch")
Correspondence = read_tsv('Data/Pairs.tsv')
Location_trees = "Trees/"

IDs = c( paste0(Correspondence$FSK, "__1"), paste0(Correspondence$FSK, "__2"), 
   paste0(Correspondence$APK, "__1"), paste0(Correspondence$APK, "__2") )
Subject = c(Correspondence$FSK, Correspondence$FSK, Correspondence$FSK, Correspondence$FSK)
Protocol = c( rep('FSK', length(Correspondence$FSK)*2 ), rep('APK', length(Correspondence$APK)*2 )  )
Repeat = c(rep(1,  length(Correspondence$FSK) ), rep(2,  length(Correspondence$FSK) ), rep(1,  length(Correspondence$FSK) ), rep(2,  length(Correspondence$FSK) )  )

Meta = tibble(ID = IDs, Subject = Subject, Protocol = Protocol, Rep = Repeat )

Tree_location = "Trees/t__SGB14546_group/RAxML_bestTree.t__SGB14546_group.StrainPhlAn4.tre"

permanova_test = function(Tree, Meta, Resample=2000, Min_N = 4, Show_plot= F, Min_to_anal=10 ){
  Meta_comm = Meta %>% filter(ID %in%  Tree$tip.label )
  KEEP = Meta_comm %>% group_by(Subject) %>% summarise(N = n()) %>% filter(N==Min_N)
  Meta_comm %>% filter(Subject %in% KEEP$Subject ) -> Meta_comm
  N = length(unique(Meta_comm$Subject))
  if (N  < Min_to_anal ){ 
    tibble(Resample = Resample ,R2 =NA, P=NA, N_samples=N, P_betadistp=NA ) -> DF 
  } else {
  Tree_filtered = keep.tip(Tree, Meta_comm$ID)
  cophenetic.phylo(Tree) -> Distances
  Distances = Distances[Meta_comm$ID, Meta_comm$ID ]
  #homogeneity of variances?
  betadisp_result <- betadisper(as.dist(Distances), Meta_comm$Protocol)
  Homogeneity = anova(betadisp_result)
  #permanova
  Resu_adonis = adonis2(Distances ~ Protocol,Meta_comm, strata = Meta_comm$Subject,permutations = Resample )
  tibble(Resample = Resample ,R2 = Resu_adonis$R2[1], P=Resu_adonis$`Pr(>F)`[1], N_samples = N, P_betadistp= Homogeneity$`Pr(>F)`[1] ) -> DF
  
  Meta_comm = Meta_comm[ match(Tree_filtered$tip.label, Meta_comm$ID ), ]
  Tree_filtered$tip.label = mutate(Meta_comm, ID2 = paste0(Subject, "_", Protocol, "_", Rep)  )$ID2
  colors = ifelse( Meta_comm$Protocol == 'APK', '#61707D', '#E85D75')
  # Plot the tree
  if (Show_plot == T){
    print('plotting')
    plot(Tree_filtered, show.tip.label = TRUE, cex = 0.4, tip.color = colors, type = "fan")
  }
  }
  return(DF)
}
Distance_test = function(Meta, Tree,SGB, Resample = 100000, Min_to_anal=10, Min_N = 4, Do_pair_wilcox=F){
  Meta_comm = Meta %>% filter(ID %in%  Tree$tip.label )
  KEEP = Meta_comm %>% group_by(Subject) %>% summarise(N = n()) %>% filter(N==Min_N)
  Meta_comm %>% filter(Subject %in% KEEP$Subject ) -> Meta_comm
  N = length(unique(Meta_comm$Subject))
  
  if (N  == 0 ){ 
    Test = tibble( SGB = SGB , Resample = Resample, Stat = NA, P =NA, N = N )
    return( list(DF = tibble() ,Test = Test , Plot = NA  ) )
  }
  
  Tree_filtered = keep.tip(Tree, Meta_comm$ID)
  cophenetic.phylo(Tree_filtered) -> Distances
    
  Distances[lower.tri(Distances)] = NA
  dist_df <- as.data.frame(as.table(Distances))
  colnames(dist_df) <- c("rowname", "colname", "distance")
  as_tibble(dist_df) %>% filter(!rowname==colname) %>% filter(! is.na(distance) ) %>% mutate(rowname = as.character(rowname), colname= as.character(colname) ) -> dist_df
  sd(dist_df$distance) -> Norm_factor
  dist_df %>% mutate(distance = distance/Norm_factor) -> dist_df
  dist_df %>% left_join(Meta_comm, by = c('rowname' = 'ID') ) %>% 
      left_join(Meta_comm, by=c('colname' = 'ID'), suffix=c("_A", "_B" ) ) -> dist_df
  dist_df %>% mutate(Species =SGB ) -> dist_df
    
  dist_df = dist_df %>% mutate(Same_subject =  ifelse(Subject_A == Subject_B, 'Same_subject', 'Different_subject') ) %>% 
      mutate(Same_protocol = ifelse(Protocol_A == Protocol_B, 'Same_protocol', 'Different_protocol')) %>% 
      mutate(Split_file = ifelse(Same_protocol == "Same_protocol" & Same_subject ==  'Same_subject', 'Same_FQ', 'Different_FQ' ) )
  dist_df = dist_df %>% mutate(Group = ifelse(Split_file == "Same_FQ", Split_file,  Same_subject ) )
  
  if (N  < Min_to_anal ){ 
    Test = tibble( SGB = SGB , Resample = Resample, Stat = NA, P =NA, N = N )
    return( list(DF =tibble() ,Test = Test , Plot = NA  ) )
  } else {
    
    #1. Make plot
    dist_df %>% ggplot(aes(x=Group, y=distance )) + geom_boxplot() + theme_bw() + scale_y_log10() -> Plot_dist #+ geom_sina(alpha=0.2)
    #2. Test
    if ( length(SGB) > 1 ){ SGB = 'all' }
    ###If paired-wilcox
    if (Do_pair_wilcox == T){
      Distances_samesubject = dist_df %>% filter(Same_subject == "Same_subject" ) %>% filter( Protocol_A == Protocol_B | Rep_B == Rep_A   ) 
      Distances_samesubject %>% select(distance, Subject_A, Group ) %>% distinct(Subject_A, Group, .keep_all = T) %>% spread(Group, distance) -> Distances_samesubject
      test = wilcox.test(Distances_samesubject$Same_FQ, Distances_samesubject$Same_subject , paired = T)
      #Distances_samesubject %>% gather(Group, Distance, 2:3 ) %>%  ggplot(aes(x=Group , y=Distance )) + geom_boxplot() + theme_bw() + scale_y_log10()
      Test_results = tibble(SGB = SGB , Resample = Resample, Stat = test$statistic, P = test$p.value , N = N  )
    }else{
      ###IF permutations
      dist_df %>%  filter(Same_subject == "Same_subject") -> Distances_samesubject
      test = coin::wilcox_test( distance ~ as.factor(Same_protocol) ,  data = Distances_samesubject, distribution = approximate(nresample = Resample))
      #independence_test( distance ~ as.factor(Same_protocol) | as.factor(Subject_A)  ,  data = Distances_samesubject,teststat = "max", distribution = approximate(nresample = Resample))
      Test_results = tibble(SGB = SGB , Resample = Resample, Stat = statistic(test), P = pvalue(test), N = N  )
    }
    return( list(DF =Distances_samesubject ,Test = Test_results, Plot = Plot_dist  ) )
    
  }

}
Make_plot_tree = function(SGB ){
  files[grepl(SGB, files)] -> Tree_loc
  Tree = read.tree(Tree_loc)
  midpoint.root(Tree) -> Tree
  Tree$tip.label = Tree$tip.label %>% str_replace('MP4_', "")  %>%  str_replace('aa', '1') %>% str_replace('ab', '2')
  
  permanova_test(Tree,Meta, 100, Min_N = 4, Show_plot = T, Min_to_anal=5 )
  
}
Prepare_DF = function(Tree_location, Meta, Test = 'permanova' ){
  set.seed(123)
  
  SGB  = basename(dirname(Tree_location))
  Tree = read.tree(Tree_location)
  midpoint.root(Tree) -> Tree
  Tree$tip.label = Tree$tip.label %>% str_replace('MP4_', "")  %>%  str_replace('aa', '1') %>% str_replace('ab', '2')
  
  #Do permanova
  if ( Test == 'permanova'){
    permanova_test(Tree,Meta, 3000, Min_N = 4) %>% mutate(SGB = SGB) %>% return()
    
  } else {
    #Do permutation distances
    Res = Distance_test(Meta, Tree, SGB, Resample = 100000, Min_to_anal=10, Min_N = 4)
    
    
    return(list(DF = Res$DF, Stats = Res$Test, Res$Plot ))
  }
}

files <- list.files(
  path = Location_trees, 
  pattern = "RAxML_bestTree\\..*\\.StrainPhlAn4\\.tre$", 
  recursive = TRUE, 
  full.names = TRUE
)

Plots_loc = "Plots_out/"

#test with permanova
Results_perm = tibble()
for (Tree_loc in files){
  print(Tree_loc)
  Res = Prepare_DF(Tree_loc, Meta, Test='permanova') %>% mutate(Tree_loc = Tree_loc) 
  Results_perm %>% rbind(Res)  -> Results_perm
}
Results_perm %>% drop_na() %>% mutate(FDR = p.adjust(P, 'fdr') ) %>% arrange(P) -> Results_perm
pdf("t__SGB5809_group_protocolTree.pdf", width = 6, height = 6, units = 'in', res = 200 )
Make_plot_tree('t__SGB5809_group')
dev.off() 
pdf("t__SGB15332_group_protocolTree.pdf", width = 6, height = 6, units = 'in', res = 200 )
Make_plot_tree('t__SGB15332_group')
dev.off() 
Make_plot_tree('t__SGB15368')

Plot_r2 = Results_perm %>% ggplot(aes(x=R2, fill=FDR<0.05 )) + geom_histogram() + theme_bw() + scale_fill_manual(values = c('grey', '#3A86FF' ) ) 
ggsave('R2_protocol.pdf',Plot_r2, width = 6, height = 6)
                                                                                                                 
write_tsv(Results_perm, 'SumamryStats_permanova.tsv')


Results = tibble()
Big_DF = tibble()
#Test with distance-based permutations
for (Tree_loc in files){
  print(Tree_loc)
  Res = Prepare_DF(Tree_loc, Meta, Test='nonparam')
  Big_DF %>% rbind(Res[[1]]) -> Big_DF
  paste0(Plots_loc, Res[[2]]$SGB, "_distance.pdf") -> NameOutPlot
  if (! is.na(Res$Stats$P) ){
    ggsave(NameOutPlot,Res[[3]], height = 4, width = 4 )
  }
  Results %>% rbind( Res[[2]] %>%  mutate(Tree_loc =  Tree_loc)  ) -> Results
}
Results %>% drop_na() %>% mutate(FDR = p.adjust(P, 'fdr') ) -> Results


ResOverall =  Big_DF %>% Compare_distances(., Resample = 1e6)
Results %>% rbind( mutate(ResOverall[[1]], Tree_loc= NA, FDR=NA )  ) -> Results

write_tsv(Results, 'SumamryStats.tsv')
ggsave(paste0(Plots_loc, "Overall_distance.pdf") ,ResOverall[[2]])
