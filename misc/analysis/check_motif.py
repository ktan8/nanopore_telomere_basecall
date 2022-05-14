#!/usr/bin/env python

import sys
import argparse



def read_fasta(fp):
	'''
	Generator code stolen from Biopython. Parses a fasta file
	and yields a single fasta line on each call
	https://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
	'''
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))


def generateMotifCombinations(sequence):
	'''
	Generate all the possible sequence combinations
	for a motif of interest after shifting the frame.
	e.g. TTAGGG, TAGGGT, AGGGTT, GGGTTA, GGTTAG, GTTAGG
	'''
	seqLength = len(sequence)
	seqCombinations = []
	for i in range(seqLength):
		newString = sequence[i:seqLength] + sequence[0:i]

		seqCombinations.append(newString)

	seqCombinations_uniq = set(seqCombinations)

	return seqCombinations_uniq

def check_motif_present(sequence, motif):
	'''
	Check if a sequence has a motif, or
	a rearrangement of it.
	'''
	motif_combinations = generateMotifCombinations(motif)

	for motif in motif_combinations:
		if motif in sequence:
			return 1

	# Return a false if no motif was found
	# across all the combinations
	return 0



def check_telomeric_motif_sequence(sequence, telomere_motif = "TTAGGG", repeat_count=3):
	'''
	Check if a sequence is possibly a telomeric
	sequence. 
	'''

	long_motif = repeat_count * telomere_motif

	if check_motif_present(sequence, long_motif):
		return 1
	else:
		return 0


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('fasta_file', metavar='fasta_file', type=str,
                    help='Fasta file to analyze')
parser.add_argument('--org', help='Organism of interest: human|arabidopsis|celegans|all', 
		    choices=['human', 'arabidopsis', 'celegans', 'all'], default='human')
args = parser.parse_args()

#fasta_file = sys.argv[1]
repeat_count = 4

# Motifs for different organisms
motif_collection = []
human_motifs = ["TTAGGG", "CCCTAA", "TTAAAA", "TTTTAA", "CTTCTT", "AAGAAG", "CCCTGG", "CCAGGG"]
arabidopsis_motifs = ["TTTAGGG", "CCCTAAA", "CCTGGG", "CCAGGG"]
celegans_motifs = ["TTAGGC", "GCCTAA"]
all_motifs = ["TTAGGG", "CCCTAA", "TTAAAA", "TTTTAA", "CTTCTT", "AAGAAG", "CCCTGG", "CCAGGG",
	      "TTTAGGG", "CCCTAAA", "TTAGGC", "GCCTAA"]

if args.org == "human":
	motif_collection = human_motifs
elif args.org == "arabidopsis":
	motif_collection = arabidopsis_motifs
elif args.org == "celegans":
	motif_collection = celegans_motifs
elif args.org == "all":
	motif_collection = all_motifs




fasta = open(fasta_file, "r")
fasta_generator = read_fasta(fasta)

header = ["readname"] + motif_collection
print("\t".join(header))
for name, sequence in fasta_generator:
	name_clean = name[1:len(name)]
	result = [check_telomeric_motif_sequence(sequence, telomere_motif=motif, repeat_count=repeat_count) for motif in motif_collection]
	result_all = [name_clean] + result
	print("\t".join(map(str, result_all)))

