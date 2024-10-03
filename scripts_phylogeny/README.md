# Information regarding running StrainPhlAn
```Do_phylogenies.py``` is the main script in this folder, and it is but the wrapper used in the HPC for running this code:
```bash
strainphlan -o {OUTPUT} -n {CPU}  --sample_with_n_markers 50 --secondary_sample_with_n_markers 50  --sample_with_n_markers_after_filt 33 --marker_in_n_samples 50 --samples {PRIMARY} --secondary_samples {SECONDARY} -c {SGB} -d /shares/CIBIO-Storage/CM/scratch/databases/metaphlansgb_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl --treeshrink --debug 
```
Selection of primary samples, which are the ones that will be used by StrainPhlAn for selecting marker genes to be used for phylogenetic reconstruction, was done from the MetaPhlAn taxonomic abundance profiles. We required a minimal depth of coverage of 2X for primary samples. Depth of coverages are calcualted using this formula:  

$$
\text{Coverage}_{\text{SGB, Sample}} = \frac{ \left( \text{AvgReadLength}_{\text{Sample}} \times \text{ReadNumber}_{\text{Sample}} \times \text{Rel.Abundance}_{\text{SGB, Sample}} \right) }{ \text{AvgGenomeLength}_{\text{SGB}} }
$$  

The code can be found in this script:
```Create_PrimarySecondary.py```

Finally, due to the fact that >20K samples were given as input for PhyloPhlAn, the standard input approach to feed the software broke. A small change in the StrainPhlAn script was done, so that a folder is specified instead, and PhyloPhlAn checks by itself if there are samples2markers files in there. This can slight modification of the script can be found in:```strainphaln_edit.py```
