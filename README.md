# SARS-CoV-2-within-host-evolution
Workflow for quantifying within-host evolution and compartmentalization of SARS-CoV-2 populations
# Summary

This workflow quantifies the temporal and spatial dynamics of intra-host single nucleotide variants (iSNVs) in longitudinally sampled populations of SARS-CoV-2. 

Scripts produce the figures associated with [Within-host evolutionary dynamics and tissue compartmentalization during acute SARS-CoV-2 infection](https://doi.org/10.1101/2022.06.21.497047), Farjo et al. (2022). 

# Requirements and setup

This project uses R (v4.1.1) and associated packages in base R.

**Additional required R packages:**

* tidyverse (v1.3.1)
* rio (v0.5.27)
* trackviewer (v1.28.1)
* gtools (v3.9.2)
* here (v1.0.1)

All files required for analysis and all script outputs are stored in "SARS_data". To begin, download this folder and create a new R project within the directory. 

# Workflow

## 1. Relationship between sample Ct and sequence quality *(Fig 1)*

1. Run "nct_vs_coverage.R"
    * Folders and files required:
	    * user_selection/nct_vs_coverage.xlsx
    * Output: 
       * **Fig 1A**
 2. Run "dilutioncontrol.R". 
    * Folders and files required:
        * dilution_output_tables/
        * dilution_controls/ct26_dilutions.csv
        * dilution_controls/ct28_dilutions.csv
        * dilution_controls/ct23_dilutions.csv
    * Output: 
        * **Fig 1B**
  3. Run "ct_vs_snps.R"
      * Folders and files required: 
	      * naive_ct.xlsx
	      * vaccinated_ct.csv
      * Output: 
	      * **Fig 1C**
	      
## 2. Within-host diversity and iSNV tracking *(Fig 2, 6, 7, S2, S3)*

1. Run "naivevariants.R"
    * Folders and files required:
	    * naive_output_tables/
    * Output:
	    * **Fig 2A**
	    * * naive_snpcounts.csv (daily iSNV counts per participant)
	    * naive_annotations.csv (iSNVs annotated by gene function)
	    * naive_variant_tables/ (iSNV frequency tracking over time, tables imported into Prism to generate **Fig 6** and **Fig S2**)
	    * naive_count_gather.RData (daily iSNV counts per participant in a ggplot-friendly format)
	    * naive_shared.RData (summary of iSNV data)
2. Run "vaccinatedvariants.R"
    * Folders and files required: 
	    * vaccinated_output_tables/
    * Output:
	    * **Fig 2B**
	    * vaccinated_snpcounts.csv (daily iSNV counts per participant)
	    * vaccinated_annotations.csv (iSNVs annotated by gene function)
	    * vaccinated_variant_tables/ (iSNV frequency tracking over time, tables imported into Prism to generate **Fig 7** and **Fig S3**)
	    * vax_count_gather.RData (daily iSNV counts per participant in ggplot-friendly format)
	    * vaccinated_shared.RData (summary of iSNV data)
3. Run "data_analysis.R"
    * Folders and files required:
	    * naive_count_gather.RData
	    * vax_count_gather.RData
	    * naive_annorations.csv
	    * vaccinated_annotations.csv
    * Output:
	    * **Fig 2C**
	    * **Fig 2D**
	    * Statistics on iSNV counts and distributions

## 3. iSNV mapping along the SARS-CoV-2 genome *(Fig 3, S1)*

1. Run "geneSNPs.R"
    * Folders and files required:
        * naive_annotations.csv
        * vaccinated_annotations.csv
     * Output:
	     * **Fig 3**
	     * **Fig S1**


## 4. Nasal samples and environmental compartmentalization *(Fig 4, S4)*

1. Run "nasal_saliva.R"
	* Folders and files required:
		* nasal_output_tables/
		* nasal_coverage_depths.RData
		* saliva_output_tables/
		* saliva_coverage_depths.RData
	* Output:
		* comparison_list.RData (R list object of tables comparing iSNV frequencies in nasal and saliva samples for each participant)
		* comparison_tables/ (tables of nasal and saliva iSNV frequencies for each participant in .csv format)
		* nasal_saliva_freqs.csv (table of all nasal and saliva iSNV frequencies, imported into Prism to generate **Fig 4A**)
2. Run "FST.R"
   * Folders and files required:
	   * comparison_list.RData
   * Output:
	   * **Fig 4B**
	   * user_fst.RData (R list object containing within- and between-environment FST values for each individual)
	   * FST_tables/ (.csv tables of within- and between-environment FST values for each individual)
	   * FST_significance.csv (significance values for comparisons of FSTs in within- vs. between-environment pairings)
3. Run "nasalvariants.R"
   * Folders and files required:
	   * nasal_output_tables/
   * Output:
	   * nasal_snpcounts.csv (daily iSNV counts per participant)
	   * nasal_variant_tables/ (iSNV frequency tracking over time, tables imported into Prism to generate **Fig S4**)

## 5. dN/dS ratios *(Fig 5)*

1. Run "dNdS.R"
	* Folders and files required:
		* naive_output_tables/
		* vaccinated_output_tables/
	* Output:
		* **Fig 5**
		* naives.RData (merged iSNV frequency and gene annotation data per participant)
		* vaccinateds.RData (merged iSNV frequency and gene annotation data per participant)

## Extras

 * all_lineages.xlsx lists Pango lineages for each sample
 * user_selection/all_user_info.xlsx contains sampling dates for all participants
 * "color_palettes.R" contains the color palettes used to generate figures
 * "coverage_depths.R" matches coverage values with samples and generates:
	 * saliva_coverage_depths.RData
	 * nasal_coverage_depths.RData
