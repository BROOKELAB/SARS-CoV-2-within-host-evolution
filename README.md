# SARS-CoV-2-within-host-evolution
Workflow for quantifying within-host evolution and compartmentalization of SARS-CoV-2 populations
# Summary

This workflow quantifies the temporal and spatial dynamics of intra-host single nucleotide variants (iSNVs) in longitudinally sampled populations of SARS-CoV-2. 

Scripts generate the figures associated with [Within-host evolutionary dynamics and tissue compartmentalization during acute SARS-CoV-2 infection](https://doi.org/10.1101/2022.06.21.497047), Farjo et al. (2022). 

# Requirements and setup

This project uses R (v4.1.1) and associated packages in base R.

**Additional required R packages:**

* tidyverse (v1.3.1)
* rio (v0.5.27)
* rlist (v0.4.6.2)
* trackViewer (v1.28.1)
* gtools (v3.9.2)
* gdata (v2.18.0.1)
* Nmisc (v0.3.7)
* data.table (v1.14.2)
* here (v1.0.1)

All files required for analysis and all script outputs are stored in `SARS_data.zip`. To begin, download this folder and create a new R project within the directory. 

# Workflow

## 1. Low-coverage sample identification

1. Run coverage.R
   * Folders and files required:
   		* `naive_depths/`
   		* `vax_depths/`
   		* `nasal_depths/`
   * Output:
   		* `naive_depth_files.RData`
     	* `naive_depth_table.csv`
     	* `vax_depth_files.RData`
     	* `vax_depth_table.csv`
     	* `nasal_depth_files.RData`
     	* `nasal_depth_table.csv`
     			
## 2. Relationship between sample Ct and sequence quality *(Fig 1, S1, S2, Table 1, Table 2)*

1. Run `nct_vs_coverage.R`
    * Folders and files required:
	    * `nct_vs_coverage.xlsx`
    * Output: 
       * **Fig 1A**
 2. Run `dilutioncontrol.R`
    * Folders and files required:
        * `dilution_output_tables/`
    * Output: 
        * **Fig 1B**
        * **Fig S1**
        * `dilution_controls/`
        * `alpha_SNPs.RData`
  3. Run `ct_vs_snps.R`
      * Folders and files required:
		* `nasal_snpcounts.RData`
    	* `naive_ct.xlsx`
		* `all_nasal_user_info.xlsx`	 
      * Output: 
	      * **Fig 1C**
	      * **Fig S2**
4. Run `sequencing_controls.R`
	* Folders and files required:
 		* `sequencing_controls/`
  		* `naive_output_tables/`
    * Output:
    	* `runcontrol_user.csv`
     	* `runcontrol_ct.csv`
      	* `runcontrol.RData`
      	* `run_cors.RData`

## 3. Error estimates *(Fig S3)*
1. Run `error_calc.R`
	* Output:
 		* **Fig S3**

## 4. Within-host diversity and iSNV tracking *(Fig 2, 5, 6, S6, S7, S8)*

1. Run `naivevariants.R`
    * Folders and files required:
	    * `naive_output_tables/`
      * `naive_depth_files.Rdata`
      * `naive_ct.xlsx`
    * Output:
	    * **Fig 2A**
	    * `naive_snpcounts.RData` (daily iSNV counts per participant)
      * `naive_intersecting.RData` (list of all shared iSNVs per participant)
      * `naive_daily_shared.RData` (daily shared iSNV counts per participant)
      * `naive_annotations.csv` (iSNVs annotated by gene function)
      * `naive_variant_tables/` (iSNV frequency tracking over time, tables imported into Prism to generate **Fig 5** and **Fig S6**)
      * `naive_shared.RData` (summary of iSNV data)
2. Run `vaccinatedvariants.R`
    * Folders and files required: 
	    * `vaccinated_output_tables/`
      * `vax_depth_files.RData`
      * `vaccinated_ct.csv`
    * Output:
	    * **Fig 2B**
	    * `vax_snpcounts.RData` (daily iSNV counts per participant)
      * `vax_intersecting.RData` (list of all shared iSNVs per participant)
      * `vax_daily_shared.RData` (daily shared iSNV counts per participant)
      * `vax_annotations.csv` (iSNVs annotated by gene function)
      * `vaccinated_variant_tables/` (iSNV frequency tracking over time, tables imported into Prism to generate **Fig 6** and **Fig S7**)
      * `vax_shared.RData` (summary of iSNV data)
3. Run `nasalvariants.R`
    * Folders and files required: 
	    * `nasal_output_tables/`
      * `nasal_depth_files.RData`
      * `all_nasal_user_info.xlsx`
    * Output:
	    * `nasal_snpcounts.RData` (daily iSNV counts per participant)
      * `nasal_intersecting.RData` (list of all shared iSNVs per participant)
      * `nasal_daily_shared.RData` (daily shared iSNV counts per participant)
      * `nasal_variant_tables/` (iSNV frequency tracking over time, tables imported into Prism to generate **Fig S8**)
      * `vax_shared.RData` (summary of iSNV data)
4. Run `data_analysis.R`
    * Folders and files required:
	    * `naive_snpcounts.RData`
	    * `vax_snpcounts.RData`
      	* `nasal_snpcounts.RData`
		* `naive_ct.xlsx`
   		* `vaccinated_ct.csv`
    	* `naive_annotations.csv`
   		 * `vaccinated_annotations.csv`
      	* `sampleIDs_saliva/`
    * Output:
	    * **Fig 2C**
	    * **Fig 2D**
	    * Statistics on iSNV counts and distributions
        
## 5. Contamination check *(Fig S9)*

1. Run `contam_check.R`
	* Folders and files required:
 		* `naive_output_tables/`
		* `vaccinated_output_tables/`
  		* `naive_depth_files.Rdata`
		* `vax_depth_files.Rdata`
		* `run_cors.RData`
  		* `all_saliva_user_info.xlsx`
	* Output:
 		* **Fig S9**

## 6. iSNV mapping along the SARS-CoV-2 genome *(Fig 3, S4, S5)*

1. Run `geneSNPs.R`
    * Folders and files required:
        * `naive_annotations.csv`
        * `vaccinated_annotations.csv`
     * Output:
	     * **Fig 3**
	     * **Fig S4**
       * **Fig S5**


## 7. Compartmentalization analysis *(Fig 4)*

1. Run `nasal_saliva.R`
	* Folders and files required:
		* `nasal_output_tables_cutsaliva/`
		* `nasal_depth_files.RData`
  		* `sampleIDs_nasalcutsaliva/`
		* `all_nasal_user_info.xlsx`
		* `saliva_output_tables/`
		* `naive_depth_files.RData`
 		* `sampleIDs_saliva/`
  		* `naive_ct.xlsx`
	* Output:
		* **Fig 4A**
		* `comparison_list.RData` (R list object of tables comparing iSNV frequencies in nasal and saliva samples for each participant)
		* `comparison_tables/` (tables of nasal and saliva iSNV frequencies for each participant in .csv format)
		* `nasal_saliva_freqs.csv` (table of all nasal and saliva iSNV frequencies, imported into Prism to generate **Fig 4A**)
2. Run `FST.R`
   * Folders and files required:
	   * `comparison_list.RData`
     * `runcontrol.RData`
   * Output:
	   * **Fig 4B,C**
	   * `user_fst.RData` (R list object containing within- and between-environment FST values for each individual)
	   * `FST_tables/` (.csv tables of within- and between-environment FST values for each individual, imported into prism to generate **Fig 4C**)
3. Run `FST_montecarlo.R`
   * Folders and files required:
     	* `user_fst.RData`
   * Output:
  		* `FST_significance.csv` (significance values derived from Monte Carlo permutation testing)
	 

## Extras
 * `misc/SARS-CoV-2_samples.xlsx` and `misc/SRA_data.xlsx` contain raw sequence information and metadata
 * `all_saliva_user_info.xlsx` and `all_nasal_user_info.xlsx` contain dates of sampling and sample barcodes for saliva and nasal samples, respectively
 * `all_lineages.xlsx` lists Pango lineages for each sample
 * `misc/432870_day7.trim.sorted.bam` contains per-nucleotide depth data for participant 432870, day 7 (referenced in text)
 * `color_palettes.R` contains the color palettes used to generate figures

