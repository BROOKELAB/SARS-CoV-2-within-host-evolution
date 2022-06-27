#!/bin/bash

#### 0. DEFINING VARIABLES ####
# gisaid metadata downloaded on June 10th, 2022
metadata=data/metadata.tsv
# reference sequence file
# (using WIV04 for consistent with GISAID alignment)
ref=data/epi_isl_402124.fasta
# name of reference sequence
ref_name=$(grep ">" $ref | sed 's/\>//g')
ref_id=$(echo $ref_name | cut -d"|" -f 2)
# corresponding GenBank data
ref_prot="data/MN996528.gff3"
# fully gapped alignment from GISAID
# downloaded on June 2nd, 2022
seqs_gapped="msa_0527/msa_0527.fasta"
mkdir figures


#############################################
########### PHYLOGENETIC ANALYSIS ###########
#############################################

#### 1. DOWNSAMPLING OF GISAID METADATA: ####
## 1A. FILTERING METADATA ##
# Filter for: 
# 1) Human hosts 
# 2) Complete 
# 3) High coverage 
# 4) Submission date complete 
# 5) Also exclude reference 
cat \
	<(head -n 1 $metadata) \
	<(awk -F'\t' '{if ($8=="Human") print $0}' $metadata | \
		awk -F'\t' '{if ($18=="True") print $0}' | \
		awk -F'\t' '{if ($19=="True") print $0}' | \
		awk -F'\t' '{if ($4 !~ "X" && $4 !~ "x") print $0}' | \
	    awk -F'\t' '{split($4,a,"-"); if (length(a)==3) print $0}') \
    > ${metadata%.tsv}_filtered_ref.tsv

# exclude ref
grep -v $ref_id ${metadata%.tsv}_filtered_ref.tsv > ${metadata%.tsv}_filtered.tsv

## 1B. DOWNSAMPLING ##
python3 scripts/downsample.py \
	--nPerMonth 100 \
	--metadata ${metadata%.tsv}_filtered.tsv

## 1C. GETTING GISAID ACCESSIONS ##
awk -F'\t' '{print $3}' ${metadata%.tsv}_filtered_100_per_month.tsv \
	> ${metadata%.tsv}_filtered_100_per_month_acc.tsv


#### 2. BUILDING TREE: ####
## 2A. ALIGN TO REFERENCE ##
# download and rename seqs manually 
# alignment without insertions 
mafft --auto --thread 4 --keeplength --addfragments \
 ${metadata%.tsv}_filtered_100_per_month.fasta $ref \
	> ${metadata%.tsv}_filtered_100_per_month_ref_aln.fasta

## 2B. BUILD TREE ##
iqtree2 -T 4 -m GTR+G4 -czb -asr -o $ref_name \
	-s ${metadata%.tsv}_filtered_100_per_month_ref_aln.fasta

## 2C: MOLECULAR CLOCK FILTER ## 
python3 scripts/clock_filter.py \
	--seqs ${metadata%.tsv}_filtered_100_per_month_ref_aln.fasta \
	--tree ${metadata%.tsv}_filtered_100_per_month_ref_aln.fasta.treefile \
	--metadata <(tail -n +2 ${metadata%.tsv}_filtered_100_per_month.tsv) \
	--metadataIDCol 2 \
	--metadataDateCol 3 \
	--root $ref_name \
	--keep $ref_name

#### 3. PLOT PHYLOGENETIC TREES FOR EACH SUBSTITUTION ####
## 3A: GET NUCLEOTIDE POSITION FOR EACH CODON ## 
rm -rf data/substitutions_ungapped_pos.csv
while read -r line; 
do 
	# get gene and codon for each substitution we are interested in
	gene=$(echo $line | cut -d',' -f1)
	codon=$(echo $line | cut -d',' -f2)
	# get starting nucleotide of gene
	gene_start=$(awk -F'\t' '{if ($3 == "gene") print $0}' $ref_prot | \
			grep gene-$gene\; | \
			cut -f 4)
	# get the \*ungapped\* position of each nucleotide in that gene
	ungapped_codon_nucs=($(expr $codon \* 3 - 3 + $gene_start) \
		$(expr $codon \* 3 - 2 + $gene_start) \
		$(expr $codon \* 3 - 1 + $gene_start))
	echo $line | tr '\n' ',' >> data/substitutions_ungapped_pos.csv
	for pos in $ungapped_codon_nucs
	do 
		echo $pos | tr '\n' ',' >> data/substitutions_ungapped_pos.csv
	done
	# add newline character
	echo "" >> data/substitutions_ungapped_pos.csv
done < data/substitutions.csv

## 3B. GET AMINO ACID IDENTITY AT NUCLEOTIDE SITES ##
# for each substitution, get nucleotides and amino acid identity at that site
# for each genome included in the tree
# used --keeplength above, so no need to convert positions
while read -r line; 
do 
	codon=$(echo $line | awk -F',' '{print $1"_"$2}')
	p1=$(echo $line | awk -F',' '{print $5}')
	p2=$(echo $line | awk -F',' '{print $6}')
	p3=$(echo $line | awk -F',' '{print $7}')
	outfile=${metadata%.tsv}_filtered_100_per_month_ref_aln_${p1}_${p2}_${p3}.tsv
	python3 scripts/isolate_fasta_region.py \
		--seqs ${metadata%.tsv}_filtered_100_per_month_ref_aln.fasta \
		--pos $p1 $p2 $p3 > $outfile
	paste \
		<(cut -f 1 $outfile) \
		<(cut -f 2- $outfile | \
			python3 scripts/translate_codon.py -) \
		> ${metadata%.tsv}_filtered_100_per_month_ref_aln_${codon}.tsv
	# create mapping between epi isl and variant
	awk -F'\t' 'FNR==NR {split($14, voc, " "); variant[$3]=voc[2]; next}{split($1, name, "|"); print $1"\t"variant[name[2]]"\t"$2}' \
		${metadata%.tsv}_filtered_100_per_month.tsv \
		${metadata%.tsv}_filtered_100_per_month_ref_aln_${codon}.tsv \
		> ${metadata%.tsv}_filtered_100_per_month_ref_aln_VOC_${codon}.tsv
done < data/substitutions_ungapped_pos.csv

## 3B. PLOT TREE ##
while read -r line; 
do 
	codon=$(echo $line | awk -F',' '{print $1"_"$2}')
	p1=$(echo $line | awk -F',' '{print $5}')
	p2=$(echo $line | awk -F',' '{print $6}')
	p3=$(echo $line | awk -F',' '{print $7}')
	outfile=${metadata%.tsv}_filtered_100_per_month_ref_aln_${p1}_${p2}_${p3}.tsv
	# finally, plot
	ptitle=$(echo $line | awk -F',' '{print $1" "$3$2 }')
	python3 scripts/plot_tree.py \
		--metadata ${metadata%.tsv}_filtered_100_per_month_ref_aln_VOC_${codon}.tsv \
		--ref $ref_name \
		--tree ${metadata%.tsv}_filtered_100_per_month_ref_aln.fasta_clockfilter.newick \
		--treeLengthMultiplier 29891 \
		--title $ptitle
done < data/substitutions_ungapped_pos.csv




#############################################
########## SUBSTITUTION HISTOGRAMS ##########
#############################################

#### 1. GAPPED SUBSTITUTION POSITIONS ####
## 1A. GAPPED POSITION MAPPING ##
# from the gapped alignment, 
# get the gapped position for each ungapped position
grep -v ">" ${seqs_gapped%.fasta}_${ref_id}.fasta | \
	awk -v FS="" '{for (i=1;i<=NF;i++) if ($i != "-") printf i"\n"}' \
	>  ${seqs_gapped%.fasta}_${ref_id}_gappedPos.txt


## 1B. GETTING NUCLEOTIDE POSITIONS FOR SUBSTITUTIONS OF INTEREST ##
rm -rf data/substitutions_gapped_pos.csv
while read -r line; 
do 
	# get gene and codon for each substitution we are interested in
	gene=$(echo $line | cut -d',' -f1)
	codon=$(echo $line | cut -d',' -f2)
	# get starting nucleotide of gene
	gene_start=$(awk -F'\t' '{if ($3 == "gene") print $0}' $ref_prot | \
			grep gene-$gene\; | \
			cut -f 4)
	# get the \*ungapped\* position of each nucleotide in that gene
	ungapped_codon_nucs=($(expr $codon \* 3 - 3 + $gene_start) \
		$(expr $codon \* 3 - 2 + $gene_start) \
		$(expr $codon \* 3 - 1 + $gene_start))
	#echo "\n" >> data/substitutions_gapped_pos.csv
	# convert each of these to a position in the gapped alignment and save
	echo $line | tr '\n' ',' >> data/substitutions_gapped_pos.csv
	for pos in $ungapped_codon_nucs
	do 
		echo $pos | tr '\n' ',' >> data/substitutions_gapped_pos.csv
	done	 
	for pos in $ungapped_codon_nucs
	do 
		sed ''"$pos"'q;d' ${seqs_gapped%.fasta}_${ref_name}_gappedPos.txt | \
			tr '\n' ',' >> data/substitutions_gapped_pos.csv
	done
	# add newline character
	echo "" >> data/substitutions_gapped_pos.csv
done < data/substitutions.csv


## 2. GET SUBSTITUTION IDENTITIES ##
# for each substitution of interest
# isolate corresponding nucleotides 
# and convert to amino acid
# for each sequence in the full alignment
# NOTE: does not account for frameshift mutations
while read -r line;
do
	codon=$(echo $line | awk -F',' '{print $1"_"$2}')
    p1=$(echo $line | awk -F',' '{print $8}')
    p2=$(echo $line | awk -F',' '{print $9}')
    p3=$(echo $line | awk -F',' '{print $10}')
    outfile=${seqs_gapped%.fasta}_${p1}_${p2}_${p3}.tsv
    python3 scripts/isolate_fasta_region.py \
    	--seqs $seqs_gapped \
    	--pos $p1 $p2 $p3  > $outfile
    paste \
		<(tail -n +2 $outfile | cut -f 1 -) \
		<(tail -n+2 $outfile | cut -f 2- -  | \
			python3 scripts/translate_codon.py -) \
		> ${outfile%.tsv}_${codon}.tsv
	# filter outfile using filtered metadata and add date
	awk -F'\t' 'FNR==NR {date[$3]=$4; next}{split($1, name, "|"); if (name[2] in date) print $1"\t"date[name[2]]"\t"$2}' \
		${metadata%.tsv}_filtered_ref.tsv \
		${outfile%.tsv}_${codon}.tsv \
		> ${outfile%.tsv}_date_${codon}.tsv 
done < data/substitutions_gapped_pos.csv

## 3. PLOT SUBSTITUTION HISTOGRAMS ##
while read -r line;
do
	codon=$(echo $line | awk -F',' '{print $1"_"$2}')
    p1=$(echo $line | awk -F',' '{print $8}')
    p2=$(echo $line | awk -F',' '{print $9}')
    p3=$(echo $line | awk -F',' '{print $10}')
    outfile=${seqs_gapped%.fasta}_${p1}_${p2}_${p3}_date_${codon}.tsv
    ptitle=$(echo $line | awk -F',' '{print $1" "$3$2 }')
	python3 scripts/plot_histogram.py \
		--metadata $outfile \
		--metadataDateCol 1 \
		--metadataAnnCol 2 \
		--ref $ref_name \
		--title $ptitle 
done < data/substitutions_gapped_pos.csv

