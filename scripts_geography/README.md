# Information regarding directory
This directory contains the scripts used for geography-phylogeny associations. This is done in the script ```Geography_vs_phylogenyMantel.R```. From there, you can obtain correlastion coefficientsh (rho values) from a mantel test, which we define as *geographical effect in the phylogeny*, which then are associated to other featues in subsequent subanalyses.
These include:  
- Association of bacterial phenotypes predicted from Traitar with *geogpraphical effect* using `GeographyEffect_vs_Features.R`. 
Prior to this Trait was ran as follows:  
```bash
traitar phenotype {location_pfam_annotation} {sample_file} from_genes {output_file}
```
PBS script to generate traitar, in addition obtain other features such as genome average genome length per SGB, and environmental prevalence of the SGB are available in ```Bacterial_Features``` 
- Script for calculating FSTs from gene alignments can be found in ```Fst/FASTA_to_FST.R```
- Asssocation between *geogpraphical effect* and strain tranmission is done in script ```Strain_sharing/4_MergeRates_correlate.R```. Scripts 1 to 3 are used to obtaine the sample-to-sample phylogenetic distances, otain the thresholds alrady available in MetaPhlAn, and generate new thresholds based on the quantiles of the distribution.

