import csv

# Store all the HGMD coding mutations into a dictionary using the Entrez gene number as the key. Record the following information:
# chromosome:position(hg19), codon number, disease name
coding_mut = {}
with open('HGMD_substitutions_2014.4.csv', 'rU') as csv_file:
	entryreader = csv.reader(csv_file)
	for entry in entryreader:
		# Entries of nonsense mutations are discarded.
		if entry[10] != "NULL" and entry[11] != "codon_number" and entry[10][-3:] not in ["TAA", "TAG", "TGA"] and entry[1] in ['DM', 'DFP']:
			# If there is no record of the gene, create a list of mutations of this gene.
			if entry[22] not in coding_mut:
				coding_mut[entry[22]] = [(entry[4], int(entry[11]), entry[8])]
			else:
				coding_mut[entry[22]].append((entry[4], int(entry[11]), entry[8]))

# Store all the mutations in a separate dictionary using a trivial key starting from 0.
mut_iter = 0
all_mut = {}
protein_num = 0
for entrez in coding_mut:
	protein_num += 1
	for mutation in coding_mut[entrez]:
		chr_pos = mutation[0]
		codon = mutation[1]
		disease = mutation[2]
		all_mut[mut_iter] = [entrez, chr_pos, codon, disease]
		mut_iter += 1

print "# of all mutations: " + str(len(all_mut))
print '# of all proteins: ' + str(protein_num)

hSIN_interaction_dict = {}
Entrez_interface_dict = {}
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
for mut_num in all_mut:
	entrez_num = all_mut[mut_num][0]
	codon_num = all_mut[mut_num][2]
	if entrez_num in valid_Entrez_list or entrez_num not in Entrez_interface_dict:
		continue
	for region in Entrez_interface_dict[entrez_num]:
		if region[0] <= codon_num <= region[1]:
			valid_Entrez_list.append(entrez_num)
			break

across_interface = []
others = []
# Find the list of mutation pairs that are across interacting interfaces. Bin them into three groups.
for mut_1 in all_mut:
	for mut_2 in all_mut:
		if mut_1 >= mut_2:
			continue
		entrez_1 = all_mut[mut_1][0]
		entrez_2 = all_mut[mut_2][0]
		codon_1 = all_mut[mut_1][2]
		codon_2 = all_mut[mut_2][2]
		if entrez_1 == entrez_2:
			continue
		if entrez_1 not in valid_Entrez_list or entrez_2 not in valid_Entrez_list:
			continue
		if (entrez_1, entrez_2) in hSIN_interaction_dict:
			interacting = False
			for interface in hSIN_interaction_dict[(entrez_1, entrez_2)]:
				if interface[0] <= codon_1 <= interface[1] and interface[2] <= codon_2 <= interface[3]:
					interacting = True
					break
			if interacting:
				across_interface.append((mut_1, mut_2))
			else:
				others.append((mut_1, mut_2))

		elif (entrez_2, entrez_1) in hSIN_interaction_dict:
			interacting = False
			for interface in hSIN_interaction_dict[(entrez_2, entrez_1)]:
				if interface[0] <= codon_2 <= interface[1] and interface[2] <= codon_1 <= interface[3]:
					interacting = True
					break
			if interacting:
				across_interface.append((mut_1, mut_2))
			else:
				others.append((mut_1, mut_2))

mapped_across_interface = 0
mapped_others = 0
same_disease_across_interface = 0
same_disease_others = 0

for pair in across_interface:
	disease_1 = all_mut[pair[0]][3].upper()
	disease_2 = all_mut[pair[1]][3].upper()
	if disease_1 not in disease_dict or disease_2 not in disease_dict:
		continue
	mapped_across_interface += 1
	if disease_dict[disease_1] == disease_dict[disease_2]:
		same_disease_across_interface += 1
for pair in others:
	disease_1 = all_mut[pair[0]][3].upper()
	disease_2 = all_mut[pair[1]][3].upper()
	if disease_1 not in disease_dict or disease_2 not in disease_dict:
		continue
	mapped_others += 1
	if disease_dict[disease_1] == disease_dict[disease_2]:
		same_disease_others += 1

print "# of mutation pairs across interfaces: " + str(mapped_across_interface)
print "# of other mutation pairs: " + str(mapped_others)
print "# of mutation pairs across interfaces causing the same disease: " + str(same_disease_across_interface)
print "# of other mutation pairs causing the same disease: " + str(same_disease_others)

