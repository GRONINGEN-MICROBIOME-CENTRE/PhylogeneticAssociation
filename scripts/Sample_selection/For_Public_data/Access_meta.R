library(tidyverse)

Merge = function(){
read_tsv("Metatada_list.txt", col_names=F) -> Metadatas

count = 0
All_df = tibble()
for (i in Metadatas$X1){
	DF = read_tsv(i)
	
	if (dim(All_df)[1] == 0){ All_df = DF 
	}else{ All_df = full_join(All_df, DF)  }
	count = count + 1
}
All_df %>% write_tsv("All_metadata.tsv")
return(All_df)
}

#All_df= Merge()
All_df = read_tsv("All_metadata.tsv")

#stool: 23522 samples
All_df %>% filter(body_site == "stool") -> DF_stool

DF_stool %>% write_tsv("Metadata_fecal.tsv")

#remove studies
#LifeLinesDeep_2016 --> I have it
Studies_out = c("LifeLinesDeep_2016", "VilaAV_2018", "WenC_2017")
#remove duplicates
read_tsv("curatedMetagenomicDataCuration/inst/extdata/duplicates_table.tsv") -> Duplicates
DF_stool %>% filter(! study_name %in% Studies_out )  %>%  filter(! (study_name == "GuptaA_2019" & sample_id  %in% Duplicates$GupDM_HAG )) -> DF_stool

#remove FMTs
DF_stool %>% filter(study_condition == "FMT") %>% group_by(study_name) %>% summarise(n()) %>% print()
Studies_out_FMT = c("IaniroG_2022 ","LiSS_2016", "DavarD_2021")
DF_stool %>% filter(! study_name %in% Studies_out_FMT ) -> DF_stool


#Remove missing
Studies_out_missingData = c("BaruchEN_2021", "DavarD_2021", "RoutyB_2018", "ContevilleLC_2019")
DF_stool %>% filter(! study_name %in% Studies_out_missingData ) -> DF_stool 


DF_stool %>% write_tsv("Metadata_fecal_selection.tsv")


