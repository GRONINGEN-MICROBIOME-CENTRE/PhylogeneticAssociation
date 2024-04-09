# Phylogenetic Association Study
This repository contain the script used for performing intraspecies associations using gut micrbobial phylogenies from dominant strains in different samples

## Selection of samples
We collected 90 publicly available studies included in https://github.com/waldronlab/curatedMetagenomicDataCuration/tree/master/inst/curated  
We added three additional Dutch publicly available studies and a Tanzanian study.
- KurilshikovA_2019 (300 Obesity, 300OB) DOI: 10.1161/CIRCRESAHA.118.314642
- GacesaR_2022 (Dutch Microbiome Project, DMP) DOI: 10.1038/s41586-022-04567-7
- ZhernakovaA_2016 (Lifelines Deep, LLD) DOI: 10.1126/science.aad3369
- StražarM_2021 (300 Tanzanian, 300TZN) DOI: 10.1038/s41467-021-25213-2

## Biobakery4 run
MetaPhLan4 and StrainPhLan4 were run with the **mpa_vJan21_CHOCOPhlAnSGB_202103.pkl** database.
The specific StrainPhlan4 script used can be found in scripts_phylogeny/  

## Associations
#### Geography vs Phylogeny
Mantel test script between genetic and geographical distances can be found in:
scripts_geography/Geography_vs_phylopgenyMantel.R
##### Strain sharing
Scripts for computing strain sharing, and its association to mantel's rho value can be bound in scripts_geography/Strain_sharing/
##### Bacterial characteristics vs Geographical effect
Scripts for associating mantel's rho value and microbial features can be found in scripts_geography/Bacterial_Features
#### Anpan associations
Scripts using for phenotype-phylogeny assocaitions (using anpan) are available in: scripts_phenotype/



 



