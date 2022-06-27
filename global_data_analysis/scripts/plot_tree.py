import baltic as bt
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from collections import defaultdict
try: 
	from utils import plot_style
except:
	from scripts.utils import plot_style

'''
THE ROOT OF THE TREE SHOULD BE AT AN X VALUE OF 0, NOT WHATEVER ARBITRARY VALUE IT IS AT
add scatter legend
'''

def plot_tree(tree_path, ax, ann_dict=None, colors_dict=None, x_attr=lambda k: k.x, 
  title=None, x_min=None):
	tree = bt.loadNewick(tree_path)
	y_dict = {i.name: i.y for i in tree.getExternal()}
	xPos = [x_attr(k) for k in tree.Objects]
	xMin = min(xPos)
	xMax = max(xPos)
	x_delta = 0
	if x_min is not None:
		x_delta = x_min - xMin
		xMin = x_min
	xRange = xMax - xMin
	tree.plotTree(ax, x_attr=lambda k: x_attr(k) + x_delta, width=0.25, colour='#333333')
	tree.plotPoints(ax, x_attr=lambda k: x_attr(k)+ x_delta, 
		colour=lambda k: colors_dict[ann_dict[k.name]],
		outline_colour='#333333', 
		size=lambda k: 100 if ann_dict[k.name] else 0, zorder=4)
	# get all used annotations
	all_aa = set([ann_dict[k.name] for k in tree.getExternal()])
	patches = [plt.plot([],[], marker="o", ms=15, markeredgewidth=2 , 
		ls="", markerfacecolor=colors_dict[aa], 
		markeredgecolor='#333333', label=aa)[0] for aa in all_aa]
	ax.legend(handles=patches, frameon=False, title=title,
		prop={'size':18}, title_fontproperties={'size':18})
	ax.set_xlim(xMin - xRange*0.025, xMax + xRange*0.025)
	ax.set_ylim(-tree.ySpan*0.025, tree.ySpan*1.025)
	ax.set_xlabel('Substitutions')
	ax.set_yticks([])
	ax.grid(axis='x', color='#eaeaea')
	_ = [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
	# add labels
	'''
	tree_annotations = sorted(set([ann_dict[k.name] for k in tree.getExternal()]) - set([None]))
	if title:
		ax.text(0.75, 0.94, title, color='#4d4d4d', size=24,
			path_effects=[path_effects.Stroke(linewidth=1, foreground='#4d4d4d'),
				path_effects.Stroke(linewidth=0.5, foreground='#333333')],
				transform=ax.transAxes, ha='right')
	x_pos = 0.775
	for label in tree_annotations:
		ax.text(x_pos, 0.94, label, color=colors_dict[label], size=24,
		path_effects=[path_effects.Stroke(linewidth=1, foreground='#333333'),
			path_effects.Stroke(linewidth=0.5, foreground=colors_dict[label])],
			transform=ax.transAxes)
		x_pos += 0.05
	'''
	return(ax, y_dict)





def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--tree', default=None)
	parser.add_argument('--colors', default='data/colors.csv')
	parser.add_argument('--metadata', default=None)
	parser.add_argument('--metadataNameCol', default=0, type=int)
	parser.add_argument('--metadataGroupCol', default=1, type=int)
	parser.add_argument('--metadataAnnCol', default=2, type=int)
	parser.add_argument('--ref')
	parser.add_argument('--title')
	parser.add_argument('--labelVOCs', nargs='+', default=['Alpha', 'Beta', 'Gamma', 'Omicron'])
	parser.add_argument('--treeLengthMultiplier', default=1, type=int)
	args = parser.parse_args()
	#args.metadata = 'data/metadata_filtered_100_per_month_ref_aln_VOC_N_389.tsv'
	#args.ref = 'hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30'
	#args.tree = 'data/metadata_filtered_100_per_month_ref_aln.fasta_clockfilter.newick'
	#args.treeLengthMultiplier = 29891
	metadata = pd.read_csv(args.metadata, sep='\t')
	ref_state = \
		metadata[metadata.iloc[:,args.metadataNameCol] == args.ref].iloc[:,args.metadataAnnCol].values[0]

	ann_dict = defaultdict(lambda: None)
	ann_dict.update({i[args.metadataNameCol]: i[args.metadataAnnCol] if 
		i[args.metadataAnnCol] != ref_state and 
			i[args.metadataAnnCol] != '-' and 
			i[args.metadataAnnCol] and 
			not pd.isna(i[args.metadataAnnCol]) else 
		None for i in metadata.values})

	colors = pd.read_csv(args.colors, sep=',', header=None)
	colors_dict=defaultdict(lambda: '#333333')
	colors_dict.update({i[0]: i[1] for i in colors.values})
	plot_style()	
	fig, ax = plt.subplots(figsize=(6.4, 4.8*2),
		constrained_layout=True)
	# setting minimum x value so root lies at 0
	ax, y_dict = plot_tree(args.tree, ax, 
		x_attr=lambda k: k.x*args.treeLengthMultiplier if k.x else 0, 
		ann_dict=ann_dict, colors_dict=colors_dict, 
		title=args.title, x_min=0)
	# for each VOC, get mean value
	label_locs = \
		np.array([(g, g_dat['seq'].map(y_dict).mean()) for 
			g, g_dat in metadata.groupby(metadata.columns[args.metadataGroupCol]) if 
			g in args.labelVOCs])
	label_locs = label_locs[label_locs[:,1].astype(float).argsort()]
	for idx, i in enumerate(label_locs):
		ax.text(ax.get_xlim()[1] + (ax.get_xlim()[1] - ax.get_xlim()[0])*0.05*(idx%2==0) , 
			i[1].astype(float), i[0], 
			rotation=-90, va='center', size=16)


	fig.savefig('figures/'+args.tree.split('/')[-1]+'_' + args.title+'.pdf')
	plt.close()
	


if __name__ == "__main__":
    run()


