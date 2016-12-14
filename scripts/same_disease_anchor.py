# Code used for the analyses in Fig 4e.

import math
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
        new_entry.append({})
        valid_mut_count += 1
        valid_mut_dict[valid_mut_count] = new_entry

# Read the TFBS information.

motif_dict = {}
with open("mut_scan.bed.txt") as motif_file:
    for row in motif_file:
        entry = row.split()
        chr = entry[0][3:]
        start = int(entry[1])
        end = int(entry[2])
        motif_num = entry[3]
        score = float(entry[4])
        if chr not in motif_dict:
            motif_dict[chr] = {(start, end): (motif_num, score)}
        else:
            if (start, end) not in motif_dict[chr]:
                motif_dict[chr][(start, end)] = (motif_num, score)

# Read the TF motif name information
TF_dict = {}
with open("TF_info.txt") as TF_file:
    for row in TF_file:
        TF_entry = row.split()
        TF_num = TF_entry[1]
        TF_name = TF_entry[2]
        TF_dict[TF_num] = TF_name

########################################################################################################################

TF_threshold = 6

# For each mutation determine which TF binding motifs it is in. Record the names of the corresponding TFs.
for mut_num in valid_mut_dict:
    chr = valid_mut_dict[mut_num][0]
    pos = valid_mut_dict[mut_num][1]
    for motif in motif_dict[chr]:
        score = motif_dict[chr][motif][1]
        if score > TF_threshold and motif[0] <= pos <= motif[1]:
            TF_num = motif_dict[chr][motif][0]
            TF_name = TF_dict[TF_num]
            if TF_name not in valid_mut_dict[mut_num][4]:
                valid_mut_dict[mut_num][4][TF_name] = motif[0]

# Only consider mutations that are located in some TF binding motifs.
mut_to_consider = {}
for mut_num in valid_mut_dict:
    if len(valid_mut_dict[mut_num][4]) > 0:
        mut_to_consider[mut_num] = valid_mut_dict[mut_num]
########################################################################################################################

# Read the mutation pairs in the same anchor.
same_anchor_list = []
with open("same_anchor_list.txt") as int_file:
    for row in int_file:
        pair = row.split()
        same_anchor_list.append((int(pair[0]) + 1, int(pair[1]) + 1)) # The mut_num in the file starts with 0, while here it starts with 1.
    int_file.close()


TF_dict = {}
# Read the TF information with the corresponding Uniprot numbers.
with open("all_TF.txt") as TF_file:
    for row in TF_file:
        TF_entry = row.split('\t')
	TF_dict[TF_entry[0]] = TF_entry[1][:-1]
    TF_file.close()

TF_interaction_list = []
# Read the HINT database to obtain protein-protein interactions.
with open("HomoSapiens_binary_hq.txt") as ppi_file:
    ppi_reader = ppi_file.readlines()
    for row in ppi_reader:
        ppi_entry = row.split('\t')
        if ppi_entry[0] != ppi_entry[1]:
            TF_interaction_list.append((ppi_entry[0], ppi_entry[1]))
    ppi_file.close()

same_motif_same_pos = []
same_motif_diff_pos = []
diff_motif_interact = []
diff_motif_non_interact = []
# Put interacting mutation pairs into three bins. Those located in the same type of motifs in interacting DNA regions, those that are not but are located in two types of motifs where the corresponding TFs interact, and those that do not fall into the above two categories.
for mut_pair in same_anchor_list:
    if mut_pair[0] in mut_to_consider and mut_pair[1] in mut_to_consider:
        same_motif_bool = False
	TFBS_pos1 = -1
	TFBS_pos2 = -2
        for TF in mut_to_consider[mut_pair[0]][4]:
            if TF in mut_to_consider[mut_pair[1]][4]:
                same_motif_bool = True
                TFBS_pos1 = mut_to_consider[mut_pair[0]][4][TF]
		TFBS_pos2 = mut_to_consider[mut_pair[1]][4][TF]
                break
        if same_motif_bool:
            if TFBS_pos1 == TFBS_pos2:
                same_motif_same_pos.append(mut_pair)
            else:
                same_motif_diff_pos.append(mut_pair)
        else:
	    interact = False
	    for TF1 in mut_to_consider[mut_pair[0]][4]:
                for TF2 in mut_to_consider[mut_pair[1]][4]:
                    if TF1 in TF_dict and TF2 in TF_dict:
                        Uniprot_TF1 = TF_dict[TF1]
			Uniprot_TF2 = TF_dict[TF2]
			if (Uniprot_TF1, Uniprot_TF2) in TF_interaction_list or (Uniprot_TF2, Uniprot_TF1) in TF_interaction_list:
                            interact = True
			    break
            if interact:
                diff_motif_interact.append(mut_pair)
            else:
                diff_motif_non_interact.append(mut_pair)

def calculate_same_disease_possibility(mutation_pair_list):
    same_disease_num = 0
    all_pair_num = 0
    possibility = 0

    for mut_pair in mutation_pair_list:
        all_pair_num += 1
        disease_1 = valid_mut_dict[mut_pair[0]][3]
        disease_2 = valid_mut_dict[mut_pair[1]][3]

        for disease in disease_1:
            if disease in disease_2:
                same_disease_num += 1

    possibility = float(same_disease_num) / all_pair_num
    SE = math.sqrt(possibility * (1 - possibility) / all_pair_num)
    print same_disease_num, all_pair_num, possibility, SE
    return [possibility, SE]

calculate_same_disease_possibility(same_motif_same_pos)
calculate_same_disease_possibility(same_motif_diff_pos)
calculate_same_disease_possibility(diff_motif_interact)
calculate_same_disease_possibility(diff_motif_non_interact)


