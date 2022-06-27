import argparse
import sys
from operator import itemgetter



def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs', 
		help='path to set of beast trees, assumed nexus format')
	parser.add_argument('--region', type=int, nargs=2, help='1 indexed')
	parser.add_argument('--pos', type=int, nargs='+', help='1 indexed')
	args = parser.parse_args()

	#args.fasta = 'data/Fig1Data_aln_excludeoutlier_Aug12_cons.fasta'
	get = []
	if args.pos:
		get.extend([i-1 for i in args.pos])
	if args.region:
		get.extend(range(args.region[0]-1, args.region[1]-1))
	
	if not get:
		# todo switch to stderr
		raise Exception('must provide at least one of "--region" or "--pos"')

	seq_len = 0
	get_str = "\t".join([str(i) for i in get])
	sys.stdout.write(f'seq\t{get_str}\n')
	with open(args.seqs, 'rt') as fp:
		getter = itemgetter(*get)
		seq = ''
		for line in fp:
			line = line.rstrip()
			# new sequence
			if line.startswith(">"):
				# if previous sequence, print output
				if seq:
					out_nucs = '\t'.join(getter(seq))
					sys.stdout.write(f'{name}\t{out_nucs}\n')
					seq_len += 1
				name = line[1:]
				seq = ''
			else:
				seq+=line

		out_nucs = '\t'.join(getter(seq))
		sys.stdout.write(f'{name}\t{out_nucs}\n')
		seq_len += 1



if __name__ == "__main__":
    run()



