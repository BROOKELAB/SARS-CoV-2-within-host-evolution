# Summary
This workflow helps to replicate the analyses presented in Figure SX of [Within-host evolutionary dynamics and tissue compartmentalization during acute SARS-CoV-2 infection](https://doi.org/10.1101/2022.06.21.497047), Farjo et al. (2022). 

The purpose of these analyses were to estimate the prevalence of specific amino acid substitutions observed within single hosts in globally sampled sequences available from the GISAID EpiCoV database (https://www.gisaid.org). The acknowledgement table for all sequences included in this analysis can be found [here](https://doi.org/10.55876/gis8.220621ms). Because of GISAID rules we cannot share the exact data files or intermediate files which were used in these analyses. However, we have provided a list of accession numbers for each of the two analyses and thus anyone with a free GISAID login should be able to download the necessary data. All analyses use epi_isl_402124 (WIV04) as a reference. The amino acids we're interested in are orf1ab 124, orf1ab 1150, orf1ab 5402, E 8, N 389.

`git status`

All of these analyses are written in Bash and Python 3 using Numpy, Pandas, Matplotlib, and Baltic. The `scripts/run_all.sh` file should run all analyses, assuming you have downloaded the necessary files. 

These scripts were primarily written by [Michael Martin](https://github.com/m-a-martin) who can be reached on GitHub or at [mmart59@emory.edu](mailto:mmart59@emory.edu). 

# Phylogenetic Analysis
This analysis downsamles the GISAID available sequences, builds a phylogenetic tree using IQtree2, and plots the resultant tree with tips labelled by amino acid identity. 

To run this analysis you will need the following starting files, all saved in `./data`: 

<ol>
  <li> `git status` `metadata.tsv`: The full GISAID metadata file. For these analyses the file was downloaded on June 10th, 2022. Happy to privately share the exact file we used to anyone with a GISAID login. </li>
  <li>`epi_isl_402124.fasta`: The epi_isl_402124 (WIV04) reference fasta file </li>
  <li>`MN996528.gff3`: The corresponding WIV04 protein annotations, downloaded from [GenBank](https://www.ncbi.nlm.nih.gov/nuccore/MN996528)</li>
  <li> </li>
</ol>

This analysis will conduct the following steps: 
<ol>
	<li> 1. Down sampling the GISAID metadata to 100 sequences per month </li>
	<ol>
		<li> A. Filter the GISAID metadata to include only complete, high coverage sequences from human hosts with complete sampling dates. This step also excludes WIV04, which will be added in later. This is conducted using Bash. </li>
		<li> B. Down sampling the filtered metadata file using the `scripts/downsample.py` Python3 script.</li>
		<li> C. Sub setting the down sampled metadata file to just the GISAID accessions with Bash. Accession file is available at `data/metadata_filtered_100_per_month_acc.tsv`. </li>
	</ol>
	<li> 2. Build ML phylogenetic tree with down sampled sequences </li>
	<ol>
		<li> 2A. Align sequences to WIV04 using MAFFT, removing any insertions relative to the reference </li>
		<li> 2B. Build tree using IQ-Tree with a GTR+G4 substitution model, collapsing near zero branches. The resultant newick file is available at ``data/metadata_filtered_100_per_month_ref_aln.fasta.treefile``. </li>
		<li> 2C. Using TreeTime through the `scripts/clock_filter.py` wrapper to remove any sequences outside of 4IQD of the expected molecular clock, rooting at WIV04 and forcing WIV04 to remain in the filtered tree. The filtered newick file is available at `data/metadata_filtered_100_per_month_ref_aln.fasta_clockfilter.newick` </li>
	</ol>
	<li> 3. Plot filtered phylogenetic tree </li>
		<ol>
			<li> 3A. Get the nucleotide position corresponding to each position in the codon of interest in WIV04 using Bash. The resultant file is available at `data/substitutions_ungapped_pos.csv`. Note: This method does not account for frameshift deletions. </li>
			<li> 3B. For each codon of interest, get the corresponding nucleotides in each sequence in the tree and translate that sequence to an amino acid.</li>
			<li> 3C. For each codon of interesting, plot the phylogenetic tree with tips annotated by amino acid identity IF they do not match the reference amino acid </li>
		</ol>

		



