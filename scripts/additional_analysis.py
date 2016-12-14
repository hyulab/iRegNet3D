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
"""
for pair in valid_pair_list:
	print mutation_dict[pair[0]][0], str(mutation_dict[pair[0]][1]), str(mutation_dict[pair[1]][1]), mutation_dict[pair[0]][2], mutation_dict[pair[1]][2]
	print '\n'
"""
"""
f = open('mut_interacting_info.txt', 'w')
for mut in mut_in_pair:
	chr = mutation_dict[mut][0]
	pos = str(mutation_dict[mut][1])
	disease = mutation_dict[mut][2]
	f.write(chr + '\t' + pos + '\t' + disease + '\n')
f.close()
print 'File writing complete!'
"""

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


"""
def measure(x):
	return fpu.up(lambda: sum((c.sup - c.inf for c in x), 0))

search_interval_dict = {}
search_interval_len = 0
with open('additional_2000.txt') as search_file:
	for line in search_file:
		entry = line.strip().split()
		search_interval = interval.interval[int(entry[1]), int(entry[2])]
		try:
			search_interval_dict[entry[0][3:]].append(search_interval)
		except KeyError:
			search_interval_dict[entry[0][3:]] = [search_interval]

for chr in search_interval_dict:
	search_in_chr = interval.interval.union(search_interval_dict[chr])
	search_interval_len += measure(search_in_chr)

def calculate_enrichment(threshold):
	interval_dict = {}
	plain_interval_dict = {}
	total_interval_len = 0
	with open('additional_scan_2000.txt') as scan_file:
		for line in scan_file:
			entry = line.strip().split()
			if float(entry[4]) < threshold:
				continue
			TF_interval = interval.interval[int(entry[1]), int(entry[2])]
			try:
				interval_dict[entry[0][3:]].append(TF_interval)
				plain_interval_dict[entry[0][3:]].append((int(entry[1]), int(entry[2])))
			except KeyError:
				interval_dict[entry[0][3:]] = [TF_interval]
				plain_interval_dict[entry[0][3:]] = [(int(entry[1]), int(entry[2]))]

	for chr in interval_dict:
		intervals_in_chr = interval.interval.union(interval_dict[chr])
		total_interval_len += measure(intervals_in_chr)

	p2 = float(total_interval_len) / search_interval_len
	
	mutation_counter = 0
	for mut in mut_in_pair:
		chr = mutation_dict[mut][0]
		pos = mutation_dict[mut][1]
		in_TFBS = False
		if chr not in plain_interval_dict:
			continue
		for motif in plain_interval_dict[chr]:
			if motif[0] <= pos <= motif[1]:
				in_TFBS = True
				break
		if in_TFBS:
			mutation_counter += 1
	
	p1 = float(mutation_counter) / len(mut_in_pair)

	enr = p1 / p2
	SE_log_enr = sqrt(1.0/mutation_counter - 1.0/len(mut_in_pair) + 1.0/total_interval_len - 1.0/search_interval_len)
	Z_score = log(enr) / SE_log_enr

	print 'Threshold = ' + str(threshold)
	print 'Enrichment: ' + str(enr)
	print 'Z score: ' + str(Z_score) + '\n'

for thr in range(6, 11):
	calculate_enrichment(thr)
"""


"""
mut_TFBS_dict = {}
threshold = 6
TF_added = []
with open('mut_scan.bed.txt') as scan_file:
	for line in scan_file:
		entry = line.strip().split()
		chr = entry[0][3:]
		start = int(entry[1])
		end = int(entry[2])
		TF = TFBS_dict[entry[3]]
		score = float(entry[4])
		if score < threshold:
			continue
		if TF not in TF_to_uni or TF_to_uni[TF] not in TF_with_int:
			continue
		if chr not in mut_by_chr:
			continue
		for mut in mut_by_chr[chr]:
			pos = mutation_dict[mut][1]
			if start <= pos <= end:
				if TF_to_uni[TF] not in TF_added:
					TF_added.append(TF_to_uni[TF])
				if mut not in mut_TFBS_dict:
					mut_TFBS_dict[mut] = [TF]
				elif TF not in mut_TFBS_dict[mut]:
					mut_TFBS_dict[mut].append(TF)
print '# TFs in the subnetwork: ' + str(len(TF_added))

subnetwork_ct = 0
for interaction in TF_TF:
	if interaction[0] in TF_added and interaction[1] in TF_added:
		subnetwork_ct += 1
print '# TF-TF interactions in the subnetwork: ' + str(subnetwork_ct)

# Now we can finally calculate which TF pairs interact and which do not.

all_TF_pair_count = 0
interacting_TF_pair = 0

for pair in valid_pair_list:
	if pair[0] not in mut_TFBS_dict or pair[1] not in mut_TFBS_dict:
		continue
	TF_list_1 = mut_TFBS_dict[pair[0]]
	TF_list_2 = mut_TFBS_dict[pair[1]]
	#all_TF_pair_count += len(TF_list_1) * len(TF_list_2)
	all_TF_pair_count += 1
	for TF_1 in TF_list_1:
		for TF_2 in TF_list_2:
			uni_1 = TF_to_uni[TF_1]
			uni_2 = TF_to_uni[TF_2]
			if (uni_1, uni_2) in TF_TF:
				#interacting_TF_pair += 1
				interacting_TF_pair += 1.0/(len(TF_list_1) * len(TF_list_2))

print all_TF_pair_count
print interacting_TF_pair
"""



"""
anchor_region_dict = {}
anchor_search_dict = {}
with open('/Users/Charlesliang/Desktop/Cornell/HY lab/For Paper/HiC_Old_Data/All_anchors.txt') as anc_file:
	for line in anc_file.readlines()[1:]:
		entry = line.strip().split('\t')
		anchor_search_dict[entry[3]] = (int(entry[1]), int(entry[2]))
		try:
			anchor_region_dict[entry[0][3:]].append((int(entry[1]), int(entry[2]), entry[3]))
		except KeyError:
			anchor_region_dict[entry[0][3:]] = [(int(entry[1]), int(entry[2]), entry[3])]

interaction_dict = {}
with open('/Users/Charlesliang/Desktop/Cornell/HY lab/For Paper/HiC_Old_Data/Target_of_all_anchors.txt') as int_file:
	for line in int_file.readlines()[1:]:
		entry = line.strip().split('\t')
		try:
			interaction_dict[entry[3]].append((int(entry[1]), int(entry[2])))
		except KeyError:
			interaction_dict[entry[3]] = [(int(entry[1]), int(entry[2]))]

mut_in_anchor_dict = {}
for mut in mutation_dict:
	chr = mutation_dict[mut][0]
	pos = mutation_dict[mut][1]
	for anchor in anchor_region_dict[chr]:
		if anchor[0] <= pos <= anchor[1]:
			mut_in_anchor_dict[mut] = anchor[2]

with open('interacting_pair.txt', 'w+') as new_file:
	for mut_1 in mut_in_anchor_dict:
		for mut_2 in mutation_dict:
			if mut_1 < mut_2:
				if mutation_dict[mut_1][0] != mutation_dict[mut_2][0]:
					continue
				pos_1 = mutation_dict[mut_1][1]
				pos_2 = mutation_dict[mut_2][1]
				anchor_1 = mut_in_anchor_dict[mut_1]
				if anchor_search_dict[anchor_1][0] <= pos_2 <= anchor_search_dict[anchor_1][1]:
					continue
				interact = False
				if anchor_1 not in interaction_dict:
					continue
				for region in interaction_dict[anchor_1]:
					if region[0] <= pos_2 <= region[1]:
						new_file.write(str(mut_1) + '\t' + str(mut_2) + '\n')
						break
"""					