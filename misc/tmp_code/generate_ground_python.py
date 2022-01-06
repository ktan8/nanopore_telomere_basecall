#!/usr/bin/env python

import sys
import numpy as np

ground_seq_file = sys.argv[1]
ground_seq_dict = dict()

f = open(ground_seq_file, "r")

for line in f:
	[readname, seq, seq_ground] = line.strip().split("\t")
	ground_seq_dict[seq] = seq_ground



# Reference sequences 'target'
ref = np.load("references.npy")

# Length of reference sequences. 
# 1:A, 2:C, 3:G, 4:T
ref_len = np.load("reference_lengths.npy")

# Normalized Nanopore data
chunk = np.load("chunks.npy")

# Decode the numbers to sequence
transdict = {0: 'N', 1: 'A', 2:'C', 3:'G', 4:'T'}
transdict_rev = {v:k for k,v in transdict.items()}

sequence_new_np = np.empty(shape = ref.shape, dtype=int)
for i in range(ref.shape[0]):
	#print("".join([transdict[x] for x in ref[i,:]]))
	seq = "".join([transdict[x] for x in ref[i,:]])
	print(seq)
	idx_last_nt = np.where((ref[i,:] == 0))[0][0]
	seq_clean = seq[0:idx_last_nt]
	print(ground_seq_dict[seq_clean])


	# Converts the seq from one motif to another
	#seq_converted = seq.replace(original, replacement)
	#print(seq_converted)
	#seq_converted_int = [transdict_rev[x] for x in seq_converted]
	#print(",".join(map(str, seq_converted_int)))
	

	#seq_int_np = np.array(seq_converted_int)
	#print(seq_int_np)
	#sequence_new_np[i,:] = seq_int_np


