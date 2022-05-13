#!/usr/bin/env python

import sys


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


fasta_file = sys.argv[1]
repeat_count = 4


motif_collection = ["TTAGGG", "CCCTAA", "TTAAAA", "TTTTAA", "CTTCTT", "AAGAAG", "CCCTGG", "CCAGGG"]

fasta = open(fasta_file, "r")
fasta_generator = read_fasta(fasta)

header = ["readname"] + motif_collection
print("\t".join(header))
for name, sequence in fasta_generator:
	name_clean = name[1:len(name)]
	result = [check_telomeric_motif_sequence(sequence, telomere_motif=motif, repeat_count=repeat_count) for motif in motif_collection]
	result_all = [name_clean] + result
	print("\t".join(map(str, result_all)))

