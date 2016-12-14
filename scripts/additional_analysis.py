import interval
from interval import fpu
from math import sqrt, log
from sets import Set
import random

# Read all the TF information upfront.

TF_to_uni = {}
with open('/Users/Charlesliang/Desktop/Research/2016/iregnet3d/static/mapping.txt') as TF_file:
	for line in TF_file.readlines()[1:]:
		entry = line.strip().split('\t')
		TF_to_uni[entry[2]] = entry[0]

TF_TF = []
TF_with_int = []
with open('/Users/Charlesliang/Desktop/Research/2016/iregnet3d/static/TF_TF_modified.txt') as TFTF_file:
	for line in TFTF_file.readlines()[1:]:
		entry = line.strip().split('\t')
		if (entry[0], entry[1]) not in TF_TF:
			TF_TF.append((entry[0], entry[1]))
		if (entry[1], entry[0]) not in TF_TF:
			TF_TF.append((entry[1], entry[0]))
		if entry[0] not in TF_with_int:
			TF_with_int.append(entry[0])
		if entry[1] not in TF_with_int:
			TF_with_int.append(entry[1])

print '# TF with interactions: ' + str(len(TF_with_int))
print '# TF-TF interactions: ' + str(len(TF_TF))
TF_TF_set = Set(TF_TF)

# Select mutations across interacting chromatin regions.

disease_dict = {}
disease_number_list = []
with open('../Noncoding_HiC/hgmd_noncoding_disease_list.txt') as dis_file:
	for line in dis_file:
		entry = line.strip().split('\t')
		disease_dict[entry[1]] = entry[0]
with open('../Noncoding_HiC/hgmd_noncoding_parsed.txt') as map_file:
	for line in map_file:
		entry = line.strip().split('\t')
		if len(entry) == 5 and entry[0] not in disease_number_list:
			disease_number_list.append(entry[0])

mutation_dict = {}
counter = 1
with open('../regulatory-2015.2.txt') as mutation_file:
	for line in mutation_file.readlines()[1:]:
		entry = line.strip().split('\t')
		if entry[6] in ['DM', 'DM?', 'DP', 'DFP']:
			if disease_dict[entry[1]] in disease_number_list:
				mutation_dict[counter] = [entry[8], int(entry[9]), entry[1], entry[2], entry[3], entry[4], entry[5], entry[7]]
				counter += 1

valid_pair_list = []
mut_in_pair = []

with open('interacting_pair.txt') as int_file:
	for line in int_file:
		entry = line.strip().split()
		valid_pair_list.append((int(entry[0]), int(entry[1])))
		if int(entry[0]) not in mut_in_pair:
			mut_in_pair.append(int(entry[0]))
		if int(entry[1]) not in mut_in_pair:
			mut_in_pair.append(int(entry[1]))

# Time to deal with TFBS data! First read TFBS-TF mapping.

TFBS_dict = {}
with open('TF_info.txt') as TF_file:
	for line in TF_file:
		entry = line.strip().split()
		TFBS_dict[entry[1]] = entry[2]

# Then read TFBS searching results. Find TFs whose binding sites mutations are located.

mut_by_chr = {}
for mut in mut_in_pair:
	chr = mutation_dict[mut][0]
	try:
		mut_by_chr[chr].append(mut)
	except KeyError:
		mut_by_chr[chr] = [mut]

motif_dict = {}
threshold = 7
iter = 0
with open('additional_scan_2000.txt') as scan_file:
	for line in scan_file:
		entry = line.strip().split()
		chr = entry[0][3:]
		start = int(entry[1])
		end = int(entry[2])
		TF = TFBS_dict[entry[3]]
		score = float(entry[4])
		if score < threshold:
			continue
		if TF not in TF_to_uni:
			continue
		if TF_to_uni[TF] not in TF_with_int:
			continue
		uni = TF_to_uni[TF]
		motif_dict[iter] = [chr, start, end, score, uni]
		iter += 1

mut_TFBS_dict = {}
for mut in mut_in_pair:
	for motif in motif_dict:
		if motif_dict[motif][0] != mutation_dict[mut][0]:
			continue
		pos = mutation_dict[mut][1]
		if abs(motif_dict[motif][1] - pos) <= 2000 and abs(motif_dict[motif][2] - pos) <= 2000:
			try:
				mut_TFBS_dict[mut].append(motif)
			except KeyError:
				mut_TFBS_dict[mut] = [motif]
print "mut_TFBS_dict has been built!"

TFBS_pair_set = Set([])
for mut_pair in valid_pair_list:
	for motif_1 in mut_TFBS_dict[mut_pair[0]]:
		for motif_2 in mut_TFBS_dict[mut_pair[1]]:
			if motif_1 < motif_2:
				tmp_set = Set([(motif_1, motif_2)])
				TFBS_pair_set.update(tmp_set)
			elif motif_1 > motif_2:
				tmp_set = Set([(motif_2, motif_1)])
				TFBS_pair_set.update(tmp_set)

print len(TFBS_pair_set)
TF_pair_counter = 0
for TFBS_pair in TFBS_pair_set:
	TF_1 = motif_dict[TFBS_pair[0]][-1]
	TF_2 = motif_dict[TFBS_pair[1]][-1]
	if (TF_1, TF_2) in TF_TF_set or (TF_2, TF_1) in TF_TF_set:
		TF_pair_counter += 1
print TF_pair_counter
print float(TF_pair_counter) / len(TFBS_pair_set)

non_interacting_pair_list = []
for mut_1 in mut_in_pair:
	for mut_2 in mut_in_pair:
		if mut_1 < mut_2 and (mut_1, mut_2) not in valid_pair_list and (mut_2, mut_1) not in valid_pair_list:
			if mutation_dict[mut_1][0] == mutation_dict[mut_2][0] and abs(mutation_dict[mut_1][1] - mutation_dict[mut_2][1]) > 100000:
				non_interacting_pair_list.append((mut_1, mut_2))
non_interacting_pair_list = random.sample(non_interacting_pair_list, 20)

ctrl_TFBS_set = Set([])
for mut_pair in non_interacting_pair_list:
	for motif_1 in mut_TFBS_dict[mut_pair[0]]:
		for motif_2 in mut_TFBS_dict[mut_pair[1]]:
			if motif_1 < motif_2:
				tmp_set = Set([(motif_1, motif_2)])
				ctrl_TFBS_set.update(tmp_set)
			elif motif_1 > motif_2:
				tmp_set = Set([(motif_2, motif_1)])
				ctrl_TFBS_set.update(tmp_set)

print len(ctrl_TFBS_set)
ctrl_TF_pair_counter = 0
for TFBS_pair in ctrl_TFBS_set:
	TF_1 = motif_dict[TFBS_pair[0]][-1]
	TF_2 = motif_dict[TFBS_pair[1]][-1]
	if (TF_1, TF_2) in TF_TF_set or (TF_2, TF_1) in TF_TF_set:
		ctrl_TF_pair_counter += 1
print ctrl_TF_pair_counter
print float(ctrl_TF_pair_counter) / len(ctrl_TFBS_set)
