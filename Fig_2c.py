import math

# Read the HGMD mutations and store disease-causing mutations in a dictionary.
noncoding_mut_count = 0
noncoding_mut_dict = {}
disease_list = []

with open("regulatory-2015.2.txt") as hgmd_file:
    hgmd_reader = hgmd_file.readlines()
    del hgmd_reader[0]
    for row in hgmd_reader:
        hgmd_entry = row.split("\t")

        if hgmd_entry[6] in ["DM", "DM?", "DFP", "DP"]:
            hgmd_entry[9] = hgmd_entry[9][:-1]
            noncoding_mut_count += 1
            noncoding_mut_dict[noncoding_mut_count] = [hgmd_entry[8], int(hgmd_entry[9]), hgmd_entry[1]]
            if hgmd_entry[1] not in disease_list:
                disease_list.append(hgmd_entry[1])

iter = 0
disease_list_dict = {}

for disease in disease_list:
    iter += 1
    disease_list_dict[disease] = str(iter)

non_dis_dict = {}
all_dis_list = []
iter = 0
with open("hgmd_noncoding_parsed.txt") as non_dis_file:
    non_dis_reader = non_dis_file.readlines()
    for row in non_dis_reader:
        non_dis_entry = row.split("\t")
        non_dis_entry[-1] = non_dis_entry[-1][:-1]
        if len(non_dis_entry) == 5:
            if non_dis_entry[4] not in all_dis_list:
                all_dis_list.append(non_dis_entry[4])
            if non_dis_entry[0] not in non_dis_dict:
                non_dis_dict[non_dis_entry[0]] = [non_dis_entry[4]]
            else:
                non_dis_dict[non_dis_entry[0]].append(non_dis_entry[4])

valid_non_mut_dict = {}
new_iter = 0
for iter in noncoding_mut_dict:
    if disease_list_dict[noncoding_mut_dict[iter][2]] in non_dis_dict:
        new_entry = noncoding_mut_dict[iter]
        new_entry.append(non_dis_dict[disease_list_dict[noncoding_mut_dict[iter][2]]])
        valid_non_mut_dict[new_iter] = new_entry
        new_iter += 1
#############################################################################################################
# Read the Hi-C chromosomal interaction dataset, and store all interacting region pairs into a dictionary. Also, we need to establish a dictionary to store all regions involved in some interaction. In addition, we need to know which set of mutations is inside some anchor region.

int_region_dict = {}
chr_interaction_dict = {}
anchor_dict = {}

with open("All_anchors.txt") as anchors_file:
    anchors_reader = anchors_file.readlines()
    del anchors_reader[0]
    for row in anchors_reader:
        anchors_entry = row.split("\t")
        anchors_entry[3] = anchors_entry[3][:-1]
        anchor_dict[anchors_entry[3]] = [anchors_entry[0][3:], int(anchors_entry[1]), int(anchors_entry[2])]

anchor_region_dict = {}

for anchor_num in anchor_dict:
    if anchor_dict[anchor_num][0] not in anchor_region_dict:
        anchor_region_dict[anchor_dict[anchor_num][0]] = [(anchor_dict[anchor_num][1], anchor_dict[anchor_num][2], anchor_num)]
    else:
        anchor_region_dict[anchor_dict[anchor_num][0]].append((anchor_dict[anchor_num][1], anchor_dict[anchor_num][2], anchor_num))


# Important fact: the interaction data only have intra-chromosomal interactions.
interaction_dict = {}
with open("Target_of_all_anchors.txt") as interaction_file:
    interaction_reader = interaction_file.readlines()
    del interaction_reader[0]
    for row in interaction_reader:
        interaction_entry = row.split("\t")

        if interaction_entry[3] not in interaction_dict:
            interaction_dict[interaction_entry[3]] = [(int(interaction_entry[1]), int(interaction_entry[2]))]
        else:
            interaction_dict[interaction_entry[3]].append((int(interaction_entry[1]), int(interaction_entry[2])))


mut_in_anchor_dict = {}
for mutation in valid_non_mut_dict:
    for anchor_region in anchor_region_dict[valid_non_mut_dict[mutation][0]]:
        if valid_non_mut_dict[mutation][1] >= anchor_region[0] and valid_non_mut_dict[mutation][1] <= anchor_region[1]:
            if anchor_region[2] in interaction_dict:
                mut_in_anchor_dict[mutation] = valid_non_mut_dict[mutation]
                mut_in_anchor_dict[mutation].append(anchor_region[2])
                break

# NOW WE LOOK AT MUTATION PAIRS WHERE AT LEAST ONE OF THE MUTATION IS IN SOME ANCHOR.
same_chr_list = []
diff_chr_list = []

for mut_1 in mut_in_anchor_dict:
    for mut_2 in valid_non_mut_dict:
        if mut_2 in mut_in_anchor_dict:
            if mut_1 < mut_2:
                if mut_in_anchor_dict[mut_1][0] == mut_in_anchor_dict[mut_2][0]:
                    same_chr_list.append((mut_1, mut_2))
                else:
                    diff_chr_list.append((mut_1, mut_2))
        else:
            if mut_in_anchor_dict[mut_1][0] == valid_non_mut_dict[mut_2][0]:
                same_chr_list.append((mut_1, mut_2))
            else:
                diff_chr_list.append((mut_1, mut_2))

# For mutation pairs on the same chromosome, check whether the DNA regions they are located at interact with each other if they are located at different anchor regions.
same_anchor_list = []
interact_list = []
non_interact_list = []

def categorize(threshold):
    same_anchor_list[:] = []
    interact_list[:] = []
    non_interact_list[:] = []

    for mut_pair in same_chr_list:
        interact = False
    
        # Case 1: both mutations are inside some anchor and are located at the same anchor region.
        if mut_pair[0] in mut_in_anchor_dict and mut_pair[1] in mut_in_anchor_dict:
            if mut_in_anchor_dict[mut_pair[0]][4] == mut_in_anchor_dict[mut_pair[1]][4]:
                same_anchor_list.append(mut_pair)
                continue
    
        # Case 2: the mutations are on interacting DNA regions.
        distance = max(valid_non_mut_dict[mut_pair[0]][1], valid_non_mut_dict[mut_pair[1]][1]) - min(valid_non_mut_dict[mut_pair[0]][1], valid_non_mut_dict[mut_pair[1]][1])
        if distance >= threshold and distance < 2000000:
            for region in interaction_dict[mut_in_anchor_dict[mut_pair[0]][4]]:
                if valid_non_mut_dict[mut_pair[1]][1] >= region[0] and valid_non_mut_dict[mut_pair[1]][1] <= region[1]:
                    interact = True
                    break
            if not interact:
                if mut_pair[1] in mut_in_anchor_dict:
                    for region in interaction_dict[mut_in_anchor_dict[mut_pair[1]][4]]:
                        if mut_in_anchor_dict[mut_pair[0]][1] >= region[0] and mut_in_anchor_dict[mut_pair[0]][1] <= region[1]:
                            interact = True
                            break
            if interact:
                interact_list.append(mut_pair)
    
        # Case 3: the mutations are on non-interacting DNA regions.
            else:
                non_interact_list.append(mut_pair)
categorize(0)

# Now we have 4 groups of mutation pairs. Calculate the possibilities that those groups of mutation pairs cause the same disease.
def calculate_same_disease_possibility(mutation_pair_list):
    same_disease_num = 0
    all_pair_num = 0
    possibility = 0
    
    for mut_pair in mutation_pair_list:
        all_pair_num += 1
        disease_1 = valid_non_mut_dict[mut_pair[0]][3]
        disease_2 = valid_non_mut_dict[mut_pair[1]][3]

        for disease in disease_1:
            if disease in disease_2:
                same_disease_num += 1

    possibility = float(same_disease_num) / all_pair_num
    SE = math.sqrt(possibility * (1 - possibility) / all_pair_num)
    print "# of mutation pairs causing the same disease: " + str(same_disease_num)
    print "# of all mutation pairs: " + str(all_pair_num)
    print "Fraction of mutation pairs causing the same disease: "+ str(possibility)
    print '\n'
    return [possibility, SE]

print "Mutation pairs in the same anchor: "
result_1 = calculate_same_disease_possibility(same_anchor_list)
print "Mutation pairs across interacting regions: "
result_2 = calculate_same_disease_possibility(interact_list)
print "Mutation pairs across non-interacting regions: "
result_3 = calculate_same_disease_possibility(non_interact_list)
print "Mutation pairs on different chromosomes: "
result_4 = calculate_same_disease_possibility(diff_chr_list)

categorize(20000)
print "Mutation pairs at least 20 kb apart across interacting regions: "
result_5 = calculate_same_disease_possibility(interact_list)
print "Mutation pairs at least 20 kb apart across non-interacting regions: "
result_6 = calculate_same_disease_possibility(non_interact_list)

categorize(50000)
print "Mutation pairs at least 50 kb apart across interacting regions: "
result_7 = calculate_same_disease_possibility(interact_list)
print "Mutation pairs at least 50 kb apart across non-interacting regions: "
result_8 = calculate_same_disease_possibility(non_interact_list)

