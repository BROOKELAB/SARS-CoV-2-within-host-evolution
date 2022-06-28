#!/usr/bin/env python3


import argparse
import sys
from operator import itemgetter



def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs', 
		help='path to set of beast trees, assumed nexus format')
	parser.add_argument('--getSeq')
	args = parser.parse_args()

	#args.fasta = 'data/Fig1Data_aln_excludeoutlier_Aug12_cons.fasta'
	
	with open(args.seqs, 'rt') as fp:
		seq = ''
		for line in fp:
			line = line.rstrip()
			# new sequence
			if line.startswith(">"):
				# if previous sequence, print output
				if seq:
					if args.getSeq in name:
						sys.stdout.write(f'>{name}\n{seq}\n')
				name = line[1:]
				seq = ''
			else:
				seq+=line
		if args.getSeq in name:
			sys.stdout.write(f'>{name}\n{seq}\n')
		


if __name__ == "__main__":
    run()



