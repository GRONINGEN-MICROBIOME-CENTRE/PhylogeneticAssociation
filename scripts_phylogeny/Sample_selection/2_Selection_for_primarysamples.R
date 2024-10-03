library(tidyverse)

print("Reading data")
Samples = read_tsv("DepthCoverage.tsv", col_names=F)
Samples2 = read_tsv("Baby_DepthCoverage.tsv", col_names=F)
Samples = rbind(Samples, Samples2)
colnames(Samples) = c("ID", "SGB", "Depth_coverage")

#Threhsold is in 2, we will give the number of samples also using other thresholds
print("Filtering Depth")
Samples %>% filter(Depth_coverage >= 2) %>% group_by(SGB) %>% summarise( N_2x = n() ) -> N2
Samples %>% filter(Depth_coverage >= 1) %>% group_by(SGB) %>% summarise( N_1x = n() ) -> N1
Samples %>% filter(Depth_coverage >= 5) %>% group_by(SGB) %>% summarise( N_5x = n() ) -> N3

print("Merging different depth thresholds")
left_join(left_join(N2, N1), N3) -> All_tests

print("Saving")
write_tsv(All_tests, "Number_samples_per_threshold.tsv")





