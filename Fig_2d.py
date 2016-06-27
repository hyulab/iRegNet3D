import numpy as np
from scipy.stats import mannwhitneyu
import sys

# Read in HGMD non-coding mutations. All the coordinates are in hg19. The entries are stored in a list. Each entry is a list containing the chromosome number of the mutation, the coordinate, and the disease it is associated with.

hgmd_mut_list = []
disease_list = []

with open("regulatory-2015.2.txt") as hgmd_file:
    hgmd_reader = hgmd_file.readlines()
    del hgmd_reader[0]
    for row in hgmd_reader:
        hgmd_entry = row.split('\t')
        # Only keep DM, DM?, DP, DFP mutations.
        if hgmd_entry[6] in ['DM', 'DM?', 'DP', 'DFP']:
            hgmd_mut_list.append([hgmd_entry[8], int(hgmd_entry[9][:-1]), hgmd_entry[1]])
            if hgmd_entry[1] not in disease_list:
                disease_list.append(hgmd_entry[1])

# Read in the parsed disease name file.

iter = 0
disease_list_dict = {}
with open("hgmd_noncoding_disease_list.txt", "w+") as file:
    for disease in disease_list:
        iter += 1
        disease_list_dict[disease] = str(iter)
        file.write(str(iter) + "\t" + disease + "\n")

noncoding_dis_dict = {}
all_dis_list = []

iter = 0
with open("hgmd_noncoding_parsed.txt") as noncoding_dis_file:
    noncoding_dis_reader = noncoding_dis_file.readlines()
    for row in noncoding_dis_reader:
        noncoding_dis_entry = row.split("\t")
        noncoding_dis_entry[-1] = noncoding_dis_entry[-1][:-1]
        if len(noncoding_dis_entry) == 5:
            if noncoding_dis_entry[4] not in all_dis_list:
                all_dis_list.append(noncoding_dis_entry[4])
            if noncoding_dis_entry[0] not in noncoding_dis_dict:
                noncoding_dis_dict[noncoding_dis_entry[0]] = [noncoding_dis_entry[4]]
            else:
                noncoding_dis_dict[noncoding_dis_entry[0]].append(noncoding_dis_entry[4])

# Only keep the mutations whose disease names are matched to some MESH numbers and append the MESH information to the entries. Store the valid non-coding mutations in a dictionary.

valid_mut_dict = {}
valid_mut_count = 0
for mutation in hgmd_mut_list:
    disease_name = mutation[2]
    if disease_list_dict[disease_name] in noncoding_dis_dict:
        new_entry = mutation
        new_entry.append(noncoding_dis_dict[disease_list_dict[disease_name]])
        valid_mut_count += 1
        valid_mut_dict[valid_mut_count] = new_entry

#############################################################################################################

# For each pair of mutations, determine if they cause the same disease.

same_disease_pair = []
diff_disease_pair = []

for num_1 in valid_mut_dict:
    for num_2 in valid_mut_dict:
        if num_1 < num_2:
            disease_list_1 = valid_mut_dict[num_1][3]
            disease_list_2 = valid_mut_dict[num_2][3]
            cause_same_disease = False
            for disease in disease_list_1:
                if disease in disease_list_2:
                    cause_same_disease = True
                    break
            if cause_same_disease:
                same_disease_pair.append((num_1, num_2))
            else:
                diff_disease_pair.append((num_1, num_2))

same_chr_same_disease = []
same_chr_diff_disease = []
diff_chr_same_disease = []
diff_chr_diff_disease = []

for pair in same_disease_pair:
    if valid_mut_dict[pair[0]][0] == valid_mut_dict[pair[1]][0]:
        same_chr_same_disease.append(pair)
    else:
        diff_chr_same_disease.append(pair)

for pair in diff_disease_pair:
    if valid_mut_dict[pair[0]][0] == valid_mut_dict[pair[1]][0]:
        same_chr_diff_disease.append(pair)
    else:
        diff_chr_diff_disease.append(pair)

#############################################################################################################

# Read in the Hi-C chromosome contact data which are in matrix form.
chr_list = ["chr" + str(i) for i in range(1, 23)]
chr_list.append("chrX")

resolution = sys.argv[1]
file_path = "./GM12878_combined/" + resolution + "_resolution_intrachromosomal/"

contact_dict = {}
all_score_list = []

for chromosome in chr_list:
    contact_dict[chromosome[3:]] = {}
    file_name = file_path + chromosome + "/MAPQGE30/" + chromosome + "_" + resolution + ".RAWobserved"
    with open(file_name) as matrix_file:
        for row in matrix_file:
            entry = row.split()
            all_score_list.append(int(entry[2][:-2]))
            if int(entry[0]) not in contact_dict[chromosome[3:]]:
                contact_dict[chromosome[3:]][int(entry[0])] = {int(entry[1]): int(entry[2][:-2])}
            else:
                contact_dict[chromosome[3:]][int(entry[0])][int(entry[1])] = int(entry[2][:-2])
"""
n, bins, patches = plt.hist(all_score_list, 50, normed=1)
histo = plt.figure()
fig.savefig(resolution + "_histo.png")
"""
#############################################################################################################

# Read the normalization vectors.
def read_norm(norm_type, resolution):
    norm_dict = {}
    
    for chromosome in chr_list:
        file_path = "./GM12878_combined/" + resolution + "_resolution_intrachromosomal/" + chromosome + "/MAPQGE30/" + chromosome + "_" + resolution + "." + norm_type + "norm"
        tmp_norm_list = []
        with open(file_path) as norm_file:
            for row in norm_file:
                tmp_norm_list.append(row[:-1])
        norm_dict[chromosome[3:]] = tmp_norm_list

    return norm_dict            

#############################################################################################################

# Compare the contact scores for mutation pairs on the same chromosome causing the same disease to scores for mutation pairs on the same chromosome causing different diseases.
resolution_dict = {"5kb":5000, "10kb":10000, "25kb":25000, "50kb":50000, "100kb":100000, "250kb":250000, "500kb":500000, "1mb":1000000}
norm_type = sys.argv[2]
norm_dict = read_norm(norm_type, resolution)

same_chr_same_disease_scores = []
for pair in same_chr_same_disease:
    chr = valid_mut_dict[pair[0]][0]
    pos_1 = valid_mut_dict[pair[0]][1]
    pos_2 = valid_mut_dict[pair[1]][1]
    if abs(pos_1 - pos_2) <= resolution_dict[resolution] or abs(pos_1 - pos_2) >= 1000000:
        continue  # Require the two mutations to be in different segments.
    if pos_1 > pos_2:
        tmp = pos_1
        pos_1 = pos_2
        pos_2 = tmp
    i_1 = pos_1 / resolution_dict[resolution]
    j_1 = pos_2 / resolution_dict[resolution]
    i = (pos_1 / resolution_dict[resolution]) * resolution_dict[resolution]
    j = (pos_2 / resolution_dict[resolution]) * resolution_dict[resolution] 

    if chr in contact_dict:
        if i in contact_dict[chr]:
            if j in contact_dict[chr][i]:
                if norm_dict[chr][i_1] != "NaN" and norm_dict[chr][j_1] != "NaN":
                    same_chr_same_disease_scores.append(contact_dict[chr][i][j] / (float(norm_dict[chr][i_1]) * float(norm_dict[chr][j_1])))
            else:
                same_chr_same_disease_scores.append(0)
        else:
            same_chr_same_disease_scores.append(0)
    else:
        same_chr_same_disease_scores.append(0)

same_chr_diff_disease_scores = []
for pair in same_chr_diff_disease:
    chr = valid_mut_dict[pair[0]][0]
    pos_1 = valid_mut_dict[pair[0]][1]
    pos_2 = valid_mut_dict[pair[1]][1]
    if abs(pos_1 - pos_2) <= resolution_dict[resolution]:
        continue
    if pos_1 > pos_2:
        tmp = pos_1
        pos_1 = pos_2
        pos_2 = tmp
    i_1 = pos_1 / resolution_dict[resolution]
    j_1 = pos_2 / resolution_dict[resolution]
    i = (pos_1 / resolution_dict[resolution]) * resolution_dict[resolution]
    j = (pos_2 / resolution_dict[resolution]) * resolution_dict[resolution]
    
    if chr in contact_dict:
        if i in contact_dict[chr]:
            if j in contact_dict[chr][i]:
                if norm_dict[chr][i_1] != "NaN" and norm_dict[chr][j_1] != "NaN":
                    same_chr_diff_disease_scores.append(contact_dict[chr][i][j] / (float(norm_dict[chr][i_1]) * float(norm_dict[chr][j_1])))
            else:
                same_chr_diff_disease_scores.append(0)
        else:
            same_chr_diff_disease_scores.append(0)
    else:
        same_chr_diff_disease_scores.append(0)
#############################################################################################################
np.asarray(same_chr_same_disease_scores)
np.asarray(same_chr_diff_disease_scores)
print "Median number of chromosomal contacts for mutation pairs causing the same disease: " + str(np.median(same_chr_same_disease_scores))
print "Median number of chromosomal contacts for mutation pairs causing different diseases: " + str(np.median(same_chr_diff_disease_scores))

U_test = mannwhitneyu(same_chr_same_disease_scores, same_chr_diff_disease_scores, alternative='two-sided')
print "P-value from Mann-Whitney U-test: " + str(U_test.pvalue)

"""
with open(resolution + "_same_chr_same_disease_scores.txt", 'w+') as same_same_file:
    for score in same_chr_same_disease_scores:
        same_same_file.write(str(score) + '\n')

with open(resolution + "_same_chr_diff_disease_scores.txt", 'w+') as same_diff_file:
    for score in same_chr_diff_disease_scores:
        same_diff_file.write(str(score) + '\n')
"""
