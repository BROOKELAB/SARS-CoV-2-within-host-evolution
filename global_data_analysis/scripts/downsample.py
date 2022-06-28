#!/usr/bin/env python3


import pandas as pd
import argparse
import numpy as np


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--metadata')
	parser.add_argument('--delim', default='\t')
	parser.add_argument('--dateCol', default='Collection date')
	parser.add_argument('--nPerMonth', type=int, default=1)
	args = parser.parse_args()
	#args.metadata = 'test.tsv'
	dat = pd.read_csv(args.metadata, sep=args.delim,  low_memory=False)
	dat[args.dateCol] = pd.to_datetime(dat[args.dateCol])
	dat['year_month'] = dat[args.dateCol].dt.to_period('M')
	n_per_month = {month: min(month_dat.shape[0], args.nPerMonth) for 
			month, month_dat in dat.groupby('year_month')}
	# now sample
	rng = np.random.default_rng(seed=111)
	rng.random
	sampled_dat = dat.groupby('year_month').apply(lambda x: 
				x.iloc[rng.choice(x.shape[0], 
					size=n_per_month[x.name], replace=False)]).reset_index(drop=True)
	sampled_dat.to_csv(f'{".".join(args.metadata.split(".")[:-1])}_{args.nPerMonth}_per_month.tsv', 
		sep='\t', index=None)





if __name__ == "__main__":
    run()

