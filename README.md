# Phylogenetic Association Study
This repository contains the scripts used for the study "**Global human gut microbiome species genetic diversity is related to geographic location and host health**." The study employs phylogenetic reconstructions of dominant gut microbial strains to investigate geographical stratification and strain-level associations with disease.

### Code Includes:
- Scripts for tree building
- Scripts for geographic associations and factors related to geographical stratification
- Scripts for tree-phenotype associations
- Downstream analysis scripts

## Sample Selection
We included 90 publicly available studies from [curatedMetagenomicDataCuration](https://github.com/waldronlab/curatedMetagenomicDataCuration/tree/master/inst/curated). Additionally, three Dutch publicly available studies and a Tanzanian study were added:
- **KurilshikovA_2019** (300 Obesity, 300OB) [DOI: 10.1161/CIRCRESAHA.118.314642](https://doi.org/10.1161/CIRCRESAHA.118.314642)
- **GacesaR_2022** (Dutch Microbiome Project, DMP) [DOI: 10.1038/s41586-022-04567-7](https://doi.org/10.1038/s41586-022-04567-7)
- **ZhernakovaA_2016** (Lifelines Deep, LLD) [DOI: 10.1126/science.aad3369](https://doi.org/10.1126/science.aad3369)
- **StražarM_2021** (300 Tanzanian, 300TZN) [DOI: 10.1038/s41467-021-25213-2](https://doi.org/10.1038/s41467-021-25213-2)

## Biobakery4 Run
MetaPhlAn4 and StrainPhlAn4 were executed using the **mpa_vJan21_CHOCOPhlAnSGB_202103.pkl** database.  
The specific StrainPhlAn4 command:  
```bash
strainphlan -o {OUTPUT} -n {CPU} --sample_with_n_markers 50 --secondary_sample_with_n_markers 50 --sample_with_n_markers_after_filt 33 --marker_in_n_samples 50 --samples {PRIMARY} --secondary_samples {SECONDARY} -c {SGB} -d /shares/CIBIO-Storage/CM/scratch/databases/metaphlansgb_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl --treeshrink --debug
```
The specific python script preparing jobs for an HPC batch run is avaialble in: 
 ```scripts_phylogeny/ ```

## Associations
#### Geography vs Phylogeny
Mantel test script between genetic and geographical distances can be found in:
scripts_geography/Geography_vs_phylopgenyMantel.R
##### Strain sharing
Scripts for computing strain sharing, and its association to mantel's rho value can be bound in ```scripts_geography/Strain_sharing/```
##### Bacterial characteristics vs Geographical effect
Scripts for associating mantel's rho value and microbial features can be found in ```scripts_geography/Bacterial_Features```
#### Anpan associations
Scripts using for phenotype-phylogeny assocaitions (using anpan) are available in: ```scripts_phenotype/```

## Data Access
To be updated!



 



