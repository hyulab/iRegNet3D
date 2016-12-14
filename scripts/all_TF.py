import csv

# Store all the HGMD coding mutations into a dictionary using the Entrez gene number as the key. Record the following information:
# chromosome:position(hg19), codon number, disease name
coding_mut = {}
with open('HGMD_substitutions_2014.4.csv', 'rU') as csv_file:
	entryreader = csv.reader(csv_file)
	for entry in entryreader:
		# Entries of nonsense mutations are discarded.
		if entry[10] != "NULL" and entry[11] != "codon_number" and entry[10][-3:] not in ["TAA", "TAG", "TGA"]:
			# If there is no record of the gene, create a list of mutations of this gene.
			if entry[22] not in coding_mut:
				coding_mut[entry[22]] = [(entry[4], int(entry[11]), entry[8])]
			else:
				coding_mut[entry[22]].append((entry[4], int(entry[11]), entry[8]))

# Read the file containing all the DNA-binding proteins.
TF_dict = {}
with open('DNABPs_Entrez.txt') as TF_file:
	for line in TF_file:
		TF_entry = line.split()
		TF_dict[TF_entry[0]] = TF_entry[1]

# Store all the mutations that are located in TF coding regions in a separate dictionary using a trivial key starting from 0.
TF_mut = {}
TF_mut_iter = 0
TF_num = 0
for entrez in coding_mut:
	if entrez in TF_dict:
		TF_num += 1
		for coding_mutation in coding_mut[entrez]:
			chr_pos = coding_mutation[0]
			codon = coding_mutation[1]
			disease = coding_mutation[2]
			TF_mut[TF_mut_iter] = [entrez, chr_pos, codon, disease]
			TF_mut_iter += 1

Entrez_interface_dict = {}
hSIN_interaction_dict = {}
Entrez_to_Uniprot = {}
# Read the hSIN interaction file. At the same time, create a mapping from Entrez number to Uniprot number.
with open("hSIN.txt") as hSIN_file:
	hSIN_reader = hSIN_file.readlines()
	del hSIN_reader[0]
	for line in hSIN_reader:
		entry = line.split()
		if entry[0] not in Entrez_to_Uniprot:
			Entrez_to_Uniprot[entry[0]] = entry[2]
			Entrez_interface_dict[entry[0]] = [(int(entry[5]), int(entry[6]))]
		else:
			if (int(entry[5]), int(entry[6])) not in Entrez_interface_dict[entry[0]]:
				Entrez_interface_dict[entry[0]].append((int(entry[5]), int(entry[6])))
		if entry[1] not in Entrez_to_Uniprot:
			Entrez_to_Uniprot[entry[1]] = entry[3]
			Entrez_interface_dict[entry[1]] = [(int(entry[8]), int(entry[9]))]
		else:
			if (int(entry[8]), int(entry[9])) not in Entrez_interface_dict[entry[1]]:
				Entrez_interface_dict[entry[1]].append((int(entry[8]), int(entry[9])))
		if entry[0] != entry[1]:
			if (entry[0], entry[1]) not in hSIN_interaction_dict:
				hSIN_interaction_dict[(entry[0], entry[1])] = [(int(entry[5]), int(entry[6]), int(entry[8]), int(entry[9]))]
			else:
				hSIN_interaction_dict[(entry[0], entry[1])].append((int(entry[5]), int(entry[6]), int(entry[8]), int(entry[9])))

disease_dict = {}
with open('hgmd2omim_map.txt') as disease_file:
	disease_reader = disease_file.readlines()
	del disease_reader[0]
	for line in disease_reader:
		disease_entry = line.split('\t')
		if disease_entry[0] not in disease_dict:
			disease_dict[disease_entry[0]] = disease_entry[1]
			
valid_Entrez_list = []
for mut_num in TF_mut:
	entrez_num = TF_mut[mut_num][0]
	codon_num = TF_mut[mut_num][2]
	if entrez_num in valid_Entrez_list or entrez_num not in Entrez_interface_dict:
		continue
	valid_Entrez_list.append(entrez_num)

all_across_interacting_TFs = []
across_non_interacting_TFs = []
# Find the list of mutation pairs that are across interacting TFs. Bin them into three groups.
for mut_1 in TF_mut:
	for mut_2 in TF_mut:
		if mut_1 >= mut_2:
			continue
		entrez_1 = TF_mut[mut_1][0]
		entrez_2 = TF_mut[mut_2][0]
		codon_1 = TF_mut[mut_1][2]
		codon_2 = TF_mut[mut_2][2]
		if entrez_1 == entrez_2:
			continue
		if entrez_1 not in valid_Entrez_list or entrez_2 not in valid_Entrez_list:
			continue
		if (entrez_1, entrez_2) in hSIN_interaction_dict:
			all_across_interacting_TFs.append((mut_1, mut_2))
		elif (entrez_2, entrez_1) in hSIN_interaction_dict:
			all_across_interacting_TFs.append((mut_1, mut_2))
		else:
			across_non_interacting_TFs.append((mut_1, mut_2))

mapped_across_all_interacting_TFs = 0
mapped_across_non_interacting_TFs = 0
same_disease_all_across_interacting_TFs = 0
same_disease_across_non_interacting_TFs = 0

for pair in all_across_interacting_TFs:
	disease_1 = TF_mut[pair[0]][3].upper()
	disease_2 = TF_mut[pair[1]][3].upper()
	if disease_1 not in disease_dict or disease_2 not in disease_dict:
		continue
	mapped_across_all_interacting_TFs += 1
	if disease_dict[disease_1] == disease_dict[disease_2]:
		same_disease_all_across_interacting_TFs += 1

for pair in across_non_interacting_TFs:
	disease_1 = TF_mut[pair[0]][3].upper()
	disease_2 = TF_mut[pair[1]][3].upper()
	if disease_1 not in disease_dict or disease_2 not in disease_dict:
		continue
	mapped_across_non_interacting_TFs += 1
	if disease_dict[disease_1] == disease_dict[disease_2]:
		same_disease_across_non_interacting_TFs += 1

print "# of mutation pairs across interacting TFs: " + str(mapped_across_all_interacting_TFs)
print "# of mutation pairs across non-interacting TFs: " + str(mapped_across_non_interacting_TFs)
print "# of mutation pairs across interacting TFs causing the same disease: " + str(same_disease_all_across_interacting_TFs)
print "# of mutation pairs across non-interacting TFs causing the same disease: " + str(same_disease_across_non_interacting_TFs)
