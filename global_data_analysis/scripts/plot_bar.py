#!/usr/bin/env python3


import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.dates as mdates
from collections import defaultdict
try: 
	from utils import plot_style
except:
	from scripts.utils import plot_style


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--colors', default='data/colors.csv')
	parser.add_argument('--metadata', default=None)
	parser.add_argument('--metadataNameCol', default=0, type=int)
	parser.add_argument('--metadataDateCol', default=1, type=int)
	parser.add_argument('--metadataAnnCol', default=2, type=int)
	parser.add_argument('--ref')
	parser.add_argument('--title')
	parser.add_argument('--maf', type=float, default=0.01)
	args = parser.parse_args()
	#args.metadata = 'msa_0527/msa_0527_1654_1655_1656_date_orf1ab_124.tsv'
	#args.ref = "hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30"
	#args.title = 'orf1ab'
	# read in colors and format
	colors = pd.read_csv(args.colors, sep=',', header=None)
	colors_dict=defaultdict(lambda: '#333333')
	colors_dict.update({i[0]: i[1] for i in colors.values})

	metadata = pd.read_csv(args.metadata, header=None, sep='\t')
	# gisaid has the country in their full alignment
	# this is not normally in gisaid sequence names
	# so we will remove
	metadata.iloc[:,args.metadataNameCol] = \
		metadata.iloc[:,args.metadataNameCol].apply(lambda k: 
			'|'.join(k.split('|')[:3]))
	# convert date object to datetime
	metadata.iloc[:,args.metadataDateCol] = \
		pd.to_datetime(metadata.iloc[:,args.metadataDateCol])
	# add year week column
	metadata['year_month'] = metadata.iloc[:,args.metadataDateCol] .dt.to_period('M')
	# remove gaps or Nans
	metadata = metadata[(metadata[2] != '-') & (~pd.isna(metadata[2]))]
	# group by year month aa, get # of each genotype
	year_month_aa_count = \
		metadata.groupby(['year_month', metadata.columns[args.metadataAnnCol]])[2].count()
	# group by year month get % of each aa
	year_month_aa_freq = \
		100*year_month_aa_count/\
			year_month_aa_count.groupby('year_month').sum()
	# fill in missing month/aa values in index
	year_month_aa_freq = year_month_aa_freq.to_frame(name='freq').reindex(
		pd.MultiIndex.from_product([year_month_aa_freq.index.get_level_values(0).unique(), 
				year_month_aa_freq.index.get_level_values(1).unique()])).fillna(0).\
		reset_index().sort_values(by='year_month')

	# keep only if max value is > args.maf
	aa_max = year_month_aa_freq.groupby(2).max()
	aa_keep = aa_max[aa_max['freq'] > args.maf].index.values
	other_rows = year_month_aa_freq[~year_month_aa_freq[2].isin(aa_keep)].index
	year_month_aa_freq.loc[other_rows,2] = 'Other'
	# now regroupby and sum
	year_month_aa_freq = year_month_aa_freq.groupby(['year_month', 2]).sum().reset_index()
	# get reference identity
	ref_aa = metadata[metadata[args.metadataNameCol] == args.ref][args.metadataAnnCol].values[0]
	# finally, alternate amino acid frequencies
	year_month_alt_aa_freq = year_month_aa_freq[year_month_aa_freq[2] != ref_aa]
	# convert tiem peirod back to datetime (plot as 1st of month)
	year_month_alt_aa_freq['year_month'] = \
		pd.to_datetime(year_month_alt_aa_freq['year_month'].astype(str))
	# define bottoms
	bottom = pd.DataFrame([year_month_alt_aa_freq['year_month'].unique()], 
		index=['year_month']).T
	bottom['bottom'] = 0.0

	# iterate over alt aa and plot
	plot_style()
	fig, ax = plt.subplots(figsize=(6.4*2, 4.8*1.5), constrained_layout=True)
	# make sure "Other" comes last
	aa_s = year_month_alt_aa_freq[args.metadataAnnCol].unique()
	if 'Other' in aa_s:
		aa_s = [i for i in aa_s if i != 'Other'] + ['Other']
	for aa in aa_s: 
		aa_dat = year_month_alt_aa_freq[year_month_alt_aa_freq[args.metadataAnnCol] == aa]
		# quite frankly, I have no idea what unit width is here
		ax.bar(aa_dat['year_month'], aa_dat['freq'], 
			width=25, color=colors_dict[aa], 
			edgecolor='#333333',
			label=aa, zorder=2)
		bottom['bottom'] +=  aa_dat['freq']

	ax.legend(frameon=False, title=args.title, 
		prop={'size':18}, title_fontproperties={'size':18})
	ax.grid(axis='y', color='#eaeaea', zorder=0)
	ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
	ax.xaxis.set_major_locator(mdates.YearLocator())
	ax.xaxis.set_minor_formatter(mdates.DateFormatter("%b"))
	ax.xaxis.set_minor_locator(mdates.MonthLocator((1,4,7,10)))
	ax.set_ylabel('Frequency (%)')
	ax.set_xlabel('Date')
	fig.savefig('figures/'+args.metadata.split('/')[-1].split('.')[0]+'.pdf')
	plt.close()


if __name__ == "__main__":
    run()



