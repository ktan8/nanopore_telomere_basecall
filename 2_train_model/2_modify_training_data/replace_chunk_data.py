#!/usr/bin/env python

# Load and replace the original chunk data
# with the ground truth


import sys
import re
import numpy as np

ground 		= sys.argv[1]
chunk 		= sys.argv[2]
references 	= sys.argv[3]
reflen		= sys.argv[4]


# ground 		= "pipeline_test/test2.chunk_groundtruth.txt"
# chunk 		= "pipeline_test/chunks.npy"
# references 	= "pipeline_test/references.npy"
# reflen		= "pipeline_test/reference_lengths.npy"






def parse_ground(ground_file):
	ground_data = open(ground_file, "r")

	ground_dict = dict()
	for line in ground_data:
		lineList 	= line.strip().split("\t")
		raw_seq 	= lineList[2]
		# Temporary fix to bug
		try:
			ground_seq 	= lineList[3]
		except:
			ground_seq 	= "NNNN"
		ground_dict[raw_seq] = ground_seq

	return ground_dict



# Reference sequences 'target'
ref = np.load(references)

# Length of reference sequences. 
# 1:A, 2:C, 3:G, 4:T
ref_len = np.load(reflen)

# Normalized Nanopore data
chunk = np.load(chunk)

# Dictionary matching the raw chunk sequence to
# the ground chunk sequence
ground_dict = parse_ground(ground)


# Decode the numbers to sequence
transdict = {0: 'N', 1: 'A', 2:'C', 3:'G', 4:'T'}
transdict_rev = {v:k for k,v in transdict.items()}


seq_ground_collections = []
seq_ground_max_len = -1

seq_with_NNN_list = [] # keep track of sequences we want to discard with NNN
seq_len_list = []

for i in range(ref.shape[0]):

		seq_raw = "".join([transdict[x] for x in ref[i,:]])
		seq_raw_replaced = seq_raw.replace("N", "") # trim off the "N" nucleotides which are indicated for padding

		# If we are unable to find the sequence in the dict,
		# we will 
		try:
			seq_ground = ground_dict[seq_raw_replaced]
		except:
			seq_ground = "NNNNN"
		seq_ground_collections.append(seq_ground)

		# Get max length for the sequences
		seq_ground_len = len(seq_ground)
		if seq_ground_len > seq_ground_max_len:
			seq_ground_max_len = seq_ground_len

		# Keep track of all ground lengths
		seq_len_list.append(seq_ground_len)

		# Keep track of sequences which have NNNN so they
		# can be discarded later
		if re.search("N", seq_ground):
			seq_with_NNN_list.append(i)


		# Check with the ground dict for the corect sequence
		# # Converts the seq from one motif to another
		# seq_converted = seq.replace(original, replacement)


		# seq_converted_int = [transdict_rev[x] for x in seq_converted]
#       #print(",".join(map(str, seq_converted_int)))
		# seq_int_np = np.array(seq_converted_int)
		#       #print(seq_int_np)
		# sequence_new_np[i,:] = seq_int_np



#print seq_ground_collections


# Generate new matrix with the data
new_shape = (ref.shape[0], seq_ground_max_len)
sequence_new_np = np.empty(shape = new_shape, dtype=int)


for i in range(len(seq_ground_collections)):
	ground_seq = seq_ground_collections[i]
	ground_seq_len = len(ground_seq)
	seq_converted_int = [transdict_rev[x] for x in ground_seq]

	seq_int_np = np.array(seq_converted_int)


	# Store the new ground sequence into the numpy
	# array
	sequence_new_np[i,0:ground_seq_len] = seq_int_np 

#print(sequence_new_np)

seq_len_list_np = np.array(seq_len_list)


# Save the raw version which has 
np.save("mod.chunks.npy", chunk) # This is same as the input chunk data
np.save("mod.references.npy", sequence_new_np)
np.save("mod.reflen.npy", seq_len_list_np)


# Generate the clean version
# Filter out those with missing sequences or "N" in the 
# ground truth sequences
clean_chunk = np.delete(chunk, seq_with_NNN_list, 0)
clean_sequence_new_np = np.delete(sequence_new_np, seq_with_NNN_list, 0)
clean_seq_len_list_np = np.delete(seq_len_list_np, seq_with_NNN_list, 0)

# Save the cleaned version
np.save("chunks.npy", clean_chunk) # This is same as the input chunk data
np.save("references.npy", clean_sequence_new_np)
np.save("reference_lengths.npy", clean_seq_len_list_np)


