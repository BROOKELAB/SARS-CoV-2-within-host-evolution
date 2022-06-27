import sys


def run():
	import sys
	from collections import defaultdict
	codon_dict=defaultdict(lambda: "NaN")
	codon_dict.update(
		{'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 
		'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 
		'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X', 
		'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W', 
		'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 
		'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 
		'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 
		'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
		'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 
		'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 
		'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 
		'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 
		'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 
		'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
		'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 
		'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 
		'---': '-'})
	for line in sys.stdin:
		# join line
		line_joined = ''.join(line.strip().split('\t')).upper()
		# chunk line
		line_chunked = [line_joined[i:i+3] for 
			i in range(0, len(line_joined),3) if 
			i+3 <=len(line_joined)]
		print("\t".join(codon_dict[i] for i in line_chunked))


if __name__ == "__main__":
    run()