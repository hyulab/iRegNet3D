# Codes for the analysis in Fig 3a.

import matplotlib as mpl
mpl.use('Agg')
from math import sqrt
from math import log
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import interval

mpl.rc('font', **{'family':'serif', 'serif':['Garamond']})
mpl.rc('text', usetex=True)
# The read_entry_file reads a file of entries with different formatting options. The key_entry_number argument is used
# to decide what the keys are in the dictionary generated. The value can be set to -1 if it is wished that entries are
# numbered from 1 on. It can also be set to a number, so the (k+1)-th value in each entry will be selected to be the
# key. After reading the entries, the first del_init_lines entries are deleted (usually the field names), and each of
# the remaining entry is split by the split_char characeter. The tail_formatting option provides a choice of formatting
# the last element in each split entry. For each entry, the tail_formatting number of characters are deleted from the
# end.
def read_entry_file (file_name, key_entry_number, del_init_lines = 1, split_char = "\t", tail_formatting = 2):
    entry_dictionary = {}
    entry_num = 0
    with open(file_name) as file:
        file_reader = file.readlines()
        for iter in range(del_init_lines):
            del file_reader[0]
        for row in file_reader:
            entry = row.split(split_char)
            if tail_formatting:
                entry[-1] = entry[-1][:-tail_formatting]
            if key_entry_number == -1:
                key = entry_num
                entry_num += 1
            elif key_entry_number < -1 or key_entry_number >= len(entry):
                return "Error: undefined key_entry_number value!"
            else:
                key = entry.pop(key_entry_number)
            entry_dictionary[key] = entry
    return entry_dictionary

# Read the non-coding mutation data.
mutation_dict = read_entry_file("regulatory-2015.2.txt", -1)

# Classify the mutations according to the type of evidence.
num_DM, num_DMq, num_DP, num_DFP, num_FP, num_other = (0,) * 6
DM_dict, DMq_dict, DP_dict, DFP_dict, FP_dict, other_dict = [{} for i in range(6)]

for mutation_entry in mutation_dict:
    if mutation_dict[mutation_entry][6] == "DM":
        DM_dict[num_DM] = mutation_dict[mutation_entry]
        num_DM += 1
    elif mutation_dict[mutation_entry][6] == "DM?":
        DMq_dict[num_DMq] = mutation_dict[mutation_entry]
        num_DMq += 1
    elif mutation_dict[mutation_entry][6] == "DP":
        DP_dict[num_DP] = mutation_dict[mutation_entry]
        num_DP += 1
    elif mutation_dict[mutation_entry][6] == "DFP":
        DFP_dict[num_DFP] = mutation_dict[mutation_entry]
        num_DFP += 1
    elif mutation_dict[mutation_entry][6] == "FP":
        FP_dict[num_FP] = mutation_dict[mutation_entry]
        num_FP += 1
    else:
        other_dict[num_other] = mutation_dict[mutation_entry]
        num_other += 1

# Calculate the total number of residues in the chromosome.
chromosome_res_num = 0
with open("hg19.chrom.sizes") as file:
    file_reader = file.readlines()
    counter = 0
    for row in file_reader:
        if counter >= 24:
            break
        counter += 1
        entry = row.split("\t")
        if entry[0] != "chrX" and entry[0] != "chrY":
            chromosome_res_num += int(entry[1])

# The cellline_calculation function counts the number of mutations in mut_dict at different genomic regions, whose
# information is provided by region_file.
def cellline_calculation(mut_dict, celltype, seg_method):
    region_file = "./" + seg_method + "/wgEncodeAwgSegmentation" + seg_method + celltype + ".bed"
    cellline_region_dict = read_entry_file(region_file, -1, 0, "\t", 1)

    # Process the promoter dictionary and enhancer dictionary, in a way that entries are classified according to
    # the chromosome.
    region_1_by_chromosome = {}
    region_2_by_chromosome = {}

    region_1_residue_num = 0
    region_2_residue_num = 0

    if seg_method == "Chromhmm":
        for iter in range(len(cellline_region_dict)):
            if cellline_region_dict[iter][3] in ["Tss", "TssF"]:
                region_1_residue_num += int(cellline_region_dict[iter][2]) - int(cellline_region_dict[iter][1]) + 1
                if cellline_region_dict[iter][0][3:] not in region_1_by_chromosome:
                    region_1_by_chromosome[cellline_region_dict[iter][0][3:]] = [(int(cellline_region_dict[iter][1]),
                                                                              int(cellline_region_dict[iter][2]))]
                else:
                    region_1_by_chromosome[cellline_region_dict[iter][0][3:]].append((int(cellline_region_dict[iter][1]),
                                                                                  int(cellline_region_dict[iter][2])))

            elif cellline_region_dict[iter][3] in ["Enh", "EnhF"]:
                region_2_residue_num += int(cellline_region_dict[iter][2]) - int(cellline_region_dict[iter][1]) + 1
                if cellline_region_dict[iter][0][3:] not in region_2_by_chromosome:
                    region_2_by_chromosome[cellline_region_dict[iter][0][3:]] = [(int(cellline_region_dict[iter][1]),
                                                                              int(cellline_region_dict[iter][2]))]
                else:
                    region_2_by_chromosome[cellline_region_dict[iter][0][3:]].append((int(cellline_region_dict[iter][1]),
                                                                                  int(cellline_region_dict[iter][2])))
    elif seg_method == "Segway":
        for iter in range(len(cellline_region_dict)):
            if cellline_region_dict[iter][3] in ["Tss", "TssF"]:
                region_1_residue_num += int(cellline_region_dict[iter][2]) - int(cellline_region_dict[iter][1]) + 1
                if cellline_region_dict[iter][0][3:] not in region_1_by_chromosome:
                    region_1_by_chromosome[cellline_region_dict[iter][0][3:]] = [(int(cellline_region_dict[iter][1]),
                                                                                  int(cellline_region_dict[iter][2]))]
                else:
                    region_1_by_chromosome[cellline_region_dict[iter][0][3:]].append(
                        (int(cellline_region_dict[iter][1]),
                         int(cellline_region_dict[iter][2])))

            elif cellline_region_dict[iter][3] in ["Enh", "EnhF", "Enh1", "Enh2", "EnhF1", "EnhF2", "EnhF3", "EnhP", "EnhPr"]:
                region_2_residue_num += int(cellline_region_dict[iter][2]) - int(cellline_region_dict[iter][1]) + 1
                if cellline_region_dict[iter][0][3:] not in region_2_by_chromosome:
                    region_2_by_chromosome[cellline_region_dict[iter][0][3:]] = [(int(cellline_region_dict[iter][1]),
                                                                                  int(cellline_region_dict[iter][2]))]
                else:
                    region_2_by_chromosome[cellline_region_dict[iter][0][3:]].append(
                        (int(cellline_region_dict[iter][1]),
                         int(cellline_region_dict[iter][2])))

    elif seg_method == "Combined":
        for iter in range(len(cellline_region_dict)):
            if cellline_region_dict[iter][3] == 'TSS':
                region_1_residue_num += int(cellline_region_dict[iter][2]) - int(cellline_region_dict[iter][1]) + 1
                if cellline_region_dict[iter][0][3:] not in region_1_by_chromosome:
                    region_1_by_chromosome[cellline_region_dict[iter][0][3:]] = [(int(cellline_region_dict[iter][1]),
                                                                          int(cellline_region_dict[iter][2]))]
                else:
                    region_1_by_chromosome[cellline_region_dict[iter][0][3:]].append((int(cellline_region_dict[iter][1]),
                                                                              int(cellline_region_dict[iter][2])))

            elif cellline_region_dict[iter][3] in ['E']:
                region_2_residue_num += int(cellline_region_dict[iter][2]) - int(cellline_region_dict[iter][1]) + 1
                if cellline_region_dict[iter][0][3:] not in region_2_by_chromosome:
                    region_2_by_chromosome[cellline_region_dict[iter][0][3:]] = [(int(cellline_region_dict[iter][1]),
                                                                          int(cellline_region_dict[iter][2]))]
                else:
                    region_2_by_chromosome[cellline_region_dict[iter][0][3:]].append((int(cellline_region_dict[iter][1]),
                                                                              int(cellline_region_dict[iter][2])))
    else:
        print "No such method!"
        return -1

    # Go through the mutations and check whether they are in the promoter region or the enhancer region.
    num_mut = 0
    num_mut_in_region_1 = 0
    num_mut_in_region_2 = 0

    for mutation in mut_dict:
        in_promoter = False
        in_enhancer = False
        num_mut += 1
        chromosome_num = mut_dict[mutation][8]
        mut_residue = int(mut_dict[mutation][9])
        if chromosome_num != "Y":
            if chromosome_num in region_1_by_chromosome:
                for pro_region in region_1_by_chromosome[chromosome_num]:
                    if pro_region[0] <= mut_residue <= pro_region[1]:
                        num_mut_in_region_1 += 1
                        in_promoter = True
                        break
            if chromosome_num in region_2_by_chromosome:
                for enh_region in region_2_by_chromosome[chromosome_num]:
                    if enh_region[0] <= mut_residue <= enh_region[1]:
                        num_mut_in_region_2 += 1
                        in_enhancer = True
                        break

    num_mut_other = num_mut - num_mut_in_region_2 - num_mut_in_region_1
    other_residue_number = chromosome_res_num - region_1_residue_num - region_2_residue_num

    # We need to return a list of the number of mutations on different genomic regions, as well as a list of the number
    # of residues of different genomic regions.
    print [[num_mut_in_region_1, num_mut_in_region_2, num_mut_other], [region_1_residue_num, region_2_residue_num,
                                                                        other_residue_number]]
    return [[num_mut_in_region_1, num_mut_in_region_2, num_mut_other], [region_1_residue_num, region_2_residue_num,
                                                                        other_residue_number]]

cellline_list = ["Gm12878", "H1hesc", "Helas3", "Hepg2", "Huvec", "K562"]
seg_method_list = ["Chromhmm", "Segway", "Combined"]
cellline_dict = {}
def calculate_enrichment (mut_dict, seg_method):
    for cellline in cellline_list:
        cellline_dict[cellline] = cellline_calculation(mut_dict, celltype=cellline, seg_method=seg_method)

    effective_mut_num = len(mut_dict)
    pro_OR_list, enh_OR_list, other_OR_list, pro_Enr_list, enh_Enr_list, other_Enr_list = [[] for i in range(6)]
    SE_OR_pro_list, SE_OR_enh_list, SE_OR_other_list, SE_Enr_pro_list, SE_Enr_enh_list, SE_Enr_other_list = [[] for i in
                                                                                                           range(6)]
    Z_list = []
    for cellline in cellline_list:
        # Calculated odds ratios and enrichment values.
        pro_p1 = cellline_dict[cellline][0][0] / float(effective_mut_num)
        enh_p1 = cellline_dict[cellline][0][1] / float(effective_mut_num)
        other_p1 = cellline_dict[cellline][0][2] / float(effective_mut_num)
        pro_p2 = cellline_dict[cellline][1][0] / float(chromosome_res_num)
        enh_p2 = cellline_dict[cellline][1][1] / float(chromosome_res_num)
        other_p2 = cellline_dict[cellline][1][2] / float(chromosome_res_num)
        pro_OR_list.append((pro_p1 / (1 - pro_p1)) / (pro_p2 / (1 - pro_p2)))
        enh_OR_list.append((enh_p1 / (1 - enh_p1)) / (enh_p2 / (1 - enh_p2)))
        other_OR_list.append((other_p1 / (1 - other_p1)) / (other_p2 / (1 - other_p2)))
        pro_Enr_list.append(pro_p1 / pro_p2)
        enh_Enr_list.append(enh_p1 / enh_p2)
        other_Enr_list.append(other_p1 / other_p2)

        # Calculate the standard errors of log odds ratios and enrichment values.
        SE_logOR_pro = sqrt((1 / float(cellline_dict[cellline][0][0])) + (
        1 / (float(effective_mut_num) - float(cellline_dict[cellline][0][0]))) + (
                              1 / float(cellline_dict[cellline][1][0])) + (
                              1 / (float(chromosome_res_num) - float(cellline_dict[cellline][1][0]))))
        SE_logOR_enh = sqrt((1 / float(cellline_dict[cellline][0][1])) + (
        1 / (float(effective_mut_num) - float(cellline_dict[cellline][0][1]))) + (
                              1 / float(cellline_dict[cellline][1][1])) + (
                              1 / (float(chromosome_res_num) - float(cellline_dict[cellline][1][1]))))
        SE_logOR_other = sqrt((1 / float(cellline_dict[cellline][0][2])) + (
        1 / (float(effective_mut_num) - float(cellline_dict[cellline][0][2]))) + (
                                1 / float(cellline_dict[cellline][1][2])) + (
                                1 / (float(chromosome_res_num) - float(cellline_dict[cellline][1][2]))))
        SE_logEnr_pro = sqrt(float((effective_mut_num - cellline_dict[cellline][0][0])) / (
        effective_mut_num * cellline_dict[cellline][0][0]) + float(
            (chromosome_res_num - cellline_dict[cellline][1][0])) / (
                            chromosome_res_num * cellline_dict[cellline][1][0]))
        SE_logEnr_enh = sqrt(float((effective_mut_num - cellline_dict[cellline][0][1])) / (
        effective_mut_num * cellline_dict[cellline][0][1]) + float(
            (chromosome_res_num - cellline_dict[cellline][1][1])) / (
                            chromosome_res_num * cellline_dict[cellline][1][1]))
        SE_logEnr_other = sqrt(float((effective_mut_num - cellline_dict[cellline][0][2])) / (
        effective_mut_num * cellline_dict[cellline][0][2]) + float(
            (chromosome_res_num - cellline_dict[cellline][1][2])) / (
                              chromosome_res_num * cellline_dict[cellline][1][2]))
        SE_OR_pro_list.append(SE_logOR_pro * (pro_p1 / (1 - pro_p1)) / (pro_p2 / (1 - pro_p2)))
        SE_OR_enh_list.append(SE_logOR_enh * (enh_p1 / (1 - enh_p1)) / (enh_p2 / (1 - enh_p2)))
        SE_OR_other_list.append(SE_logOR_other * (other_p1 / (1 - other_p1)) / (other_p2 / (1 - other_p2)))
        SE_Enr_pro_list.append(SE_logEnr_pro * pro_p1 / pro_p2)
        SE_Enr_enh_list.append(SE_logEnr_enh * enh_p1 / enh_p2)
        SE_Enr_other_list.append(SE_logEnr_other * other_p1 / other_p2)

        Z_pro = log(pro_p1 / pro_p2) / SE_logEnr_pro
        Z_enh = log(enh_p1 / enh_p2) / SE_logEnr_enh
        Z_other = log(other_p1 / other_p2) / SE_logEnr_other
        Z_list.append([Z_pro, Z_enh, Z_other])

    # Use matplotlib to plot the enrichment with error bars of the mutations in promoters, enhancers and other regions.
    N = len(cellline_list)
    ind = np.arange(N)
    width = 0.25
    fig, ax = plt.subplots()
    pro_plot = [x - 1 for x in pro_Enr_list]
    enh_plot = [x - 1 for x in enh_Enr_list]
    other_plot = [x - 1 for x in other_Enr_list]
    rects1 = ax.bar(ind + 0.1, pro_plot, width, color='#fbb4ae', yerr=SE_Enr_pro_list)
    rects2 = ax.bar(ind + width + 0.1, enh_plot, width, color='#b3cde3', yerr=SE_Enr_enh_list)
    rects3 = ax.bar(ind + 2 * width + 0.1, other_plot, width, color='#ccebc5', yerr=SE_Enr_other_list)

    ax.axhline(0, color='black', lw=1)
    plt.gca().yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + 1))
    print Z_list
    # Mark the statistically significant ones with astericks.
    for i in range(len(cellline_list)):
        for j in range(3):
            Z_score = Z_list[i][j]
            if abs(Z_score) > 1.96:
                if abs(Z_score) > 2.577:
                    if abs(Z_score) > 3.3:
                        mark = '***'
                    else:
                        mark = '**'
                else:
                    mark = '*'
            else:
                mark = ''
            if j == 0: # Promoter
                ax.annotate(mark, xy=(0.1 + ind[i] + 0.5 * width, pro_plot[i] + SE_Enr_pro_list[i] + 0.03),
                            zorder=10, ha='center')
            elif j == 1: # Enhancer
                ax.annotate(mark, xy=(0.1 + ind[i] + 1.5 * width, enh_plot[i] + SE_Enr_enh_list[i] + 0.03),
                            zorder=10, ha='center')
            else: # Other
                ax.annotate(mark, xy=(0.1 + ind[i] + 2.5 * width, other_plot[i] + SE_Enr_other_list[i] + 0.03),
                            zorder=10, ha='center')

    loc = "upper center"
    if seg_method == "Combined":
        plt.ylim((-0.5, 2))

    ax.legend((rects1[0], rects2[0], rects3[0]), ('Promoter', 'Enhancer', 'Other'), loc=loc)
    ax.set_ylabel('Enrichment of HGMD non-coding mutations')
    ax.set_xticks(ind + 1.5 * width + 0.1)
    ax.set_xticklabels(cellline_list)
    plt.savefig(seg_method + "_new.png")

def measure(x):
    from interval import fpu
    return fpu.up(lambda: sum((c.sup - c.inf for c in x), 0))

def overall_calculation(mut_dict, seg_method):

    region_dict_list = []
    region_1_by_chromosome = {}
    region_2_by_chromosome = {}

    for cellline in cellline_list:
        region_file = "./" + seg_method + "/wgEncodeAwgSegmentation" + seg_method + cellline + ".bed"
        cellline_region_dict = read_entry_file(region_file, -1, 0, "\t", 1)
        region_dict_list.append(cellline_region_dict)

    if seg_method == "Chromhmm":
        for region_dict in region_dict_list:
            for iter in range(len(region_dict)):
                if region_dict[iter][3] in ["Tss", "TssF"]:
                    if region_dict[iter][0][3:] not in region_1_by_chromosome:
                        region_1_by_chromosome[region_dict[iter][0][3:]] = [
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2]))]
                    else:
                        region_1_by_chromosome[region_dict[iter][0][3:]].append(
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2])))

                elif region_dict[iter][3] in ["Enh", "EnhF"]:
                    if region_dict[iter][0][3:] not in region_2_by_chromosome:
                        region_2_by_chromosome[region_dict[iter][0][3:]] = [
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2]))]
                    else:
                        region_2_by_chromosome[region_dict[iter][0][3:]].append(
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2])))

    elif seg_method == "Segway":
        for region_dict in region_dict_list:
            for iter in range(len(region_dict)):
                if region_dict[iter][3] in ["Tss", "TssF"]:
                    if region_dict[iter][0][3:] not in region_1_by_chromosome:
                        region_1_by_chromosome[region_dict[iter][0][3:]] = [
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2]))]
                    else:
                        region_1_by_chromosome[region_dict[iter][0][3:]].append(
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2])))

                elif region_dict[iter][3] in ["Enh", "EnhF", "Enh1", "Enh2", "EnhF1", "EnhF2", "EnhF3", "EnhP", "EnhPr"]:
                    if region_dict[iter][0][3:] not in region_2_by_chromosome:
                        region_2_by_chromosome[region_dict[iter][0][3:]] = [
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2]))]
                    else:
                        region_2_by_chromosome[region_dict[iter][0][3:]].append(
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2])))

    elif seg_method == "Combined":
        for region_dict in region_dict_list:
            for iter in range(len(region_dict)):
                if region_dict[iter][3] in ["TSS"]:
                    if region_dict[iter][0][3:] not in region_1_by_chromosome:
                        region_1_by_chromosome[region_dict[iter][0][3:]] = [
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2]))]
                    else:
                        region_1_by_chromosome[region_dict[iter][0][3:]].append(
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2])))

                elif region_dict[iter][3] in ["E"]:
                    if region_dict[iter][0][3:] not in region_2_by_chromosome:
                        region_2_by_chromosome[region_dict[iter][0][3:]] = [
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2]))]
                    else:
                        region_2_by_chromosome[region_dict[iter][0][3:]].append(
                            (int(region_dict[iter][1]),
                             int(region_dict[iter][2])))

    else:
        print "No such method!"
        return -1

    # Go through the mutations and check whether they are in the promoter region or the enhancer region.
    num_mut = 0
    num_mut_in_region_1 = 0
    num_mut_in_region_2 = 0
    num_mut_other = 0

    for mutation in mut_dict:
        in_promoter = False
        in_enhancer = False
        num_mut += 1
        chromosome_num = mut_dict[mutation][8]
        mut_residue = int(mut_dict[mutation][9])
        if chromosome_num != "Y":
            if chromosome_num in region_1_by_chromosome:
                for pro_region in region_1_by_chromosome[chromosome_num]:
                    if pro_region[0] <= mut_residue <= pro_region[1]:
                        num_mut_in_region_1 += 1
                        in_promoter = True
                        break
            if chromosome_num in region_2_by_chromosome:
                for enh_region in region_2_by_chromosome[chromosome_num]:
                    if enh_region[0] <= mut_residue <= enh_region[1]:
                        num_mut_in_region_2 += 1
                        in_enhancer = True
                        break
        if not in_promoter and not in_enhancer:
            num_mut_other += 1

    region_1_residue_num = 0
    region_2_residue_num = 0
    region_1_and_2_residue_num = 0

    for chr in region_1_by_chromosome:
        tmp_list = []
        for region in region_1_by_chromosome[chr]:
            region_1_interval = interval.interval[region[0], region[1]]
            tmp_list.append(region_1_interval)
        tmp_interval_chr = interval.interval.union(tmp_list)
        region_1_residue_num += measure(tmp_interval_chr)

    for chr in region_2_by_chromosome:
        tmp_list = []
        tmp_2_list = []
        for region in region_2_by_chromosome[chr]:
            region_2_interval = interval.interval[region[0], region[1]]
            tmp_list.append(region_2_interval)
            tmp_2_list.append(region_2_interval)
        for region in region_1_by_chromosome[chr]:
            region_1_interval = interval.interval[region[0], region[1]]
            tmp_2_list.append(region_1_interval)
        tmp_interval_chr = interval.interval.union(tmp_list)
        tmp_2_interval_chr = interval.interval.union(tmp_2_list)
        region_2_residue_num += measure(tmp_interval_chr)
        region_1_and_2_residue_num += measure(tmp_2_interval_chr)

    other_residue_number = chromosome_res_num - region_1_and_2_residue_num

    return [[num_mut_in_region_1, num_mut_in_region_2, num_mut_other], [region_1_residue_num, region_2_residue_num,
                                                                        other_residue_number]]
seg_method_dict = {}
def calculate_overall_enrichment(mut_dict):
    for seg_method in seg_method_list:
        seg_method_dict[seg_method] = overall_calculation(mut_dict, seg_method)

    effective_mut_num = len(mut_dict)
    pro_OR_list, enh_OR_list, other_OR_list, pro_Enr_list, enh_Enr_list, other_Enr_list = [[] for i in range(6)]
    SE_OR_pro_list, SE_OR_enh_list, SE_OR_other_list, SE_Enr_pro_list, SE_Enr_enh_list, SE_Enr_other_list = [[] for i in
                                                                                                             range(6)]
    Z_list = []
    for seg_method in seg_method_list:
        # Calculated odds ratios and enrichment values.
        pro_p1 = seg_method_dict[seg_method][0][0] / float(effective_mut_num)
        enh_p1 = seg_method_dict[seg_method][0][1] / float(effective_mut_num)
        other_p1 = seg_method_dict[seg_method][0][2] / float(effective_mut_num)
        pro_p2 = seg_method_dict[seg_method][1][0] / float(chromosome_res_num)
        enh_p2 = seg_method_dict[seg_method][1][1] / float(chromosome_res_num)
        other_p2 = seg_method_dict[seg_method][1][2] / float(chromosome_res_num)
        pro_OR_list.append((pro_p1 / (1 - pro_p1)) / (pro_p2 / (1 - pro_p2)))
        enh_OR_list.append((enh_p1 / (1 - enh_p1)) / (enh_p2 / (1 - enh_p2)))
        other_OR_list.append((other_p1 / (1 - other_p1)) / (other_p2 / (1 - other_p2)))
        pro_Enr_list.append(pro_p1 / pro_p2)
        enh_Enr_list.append(enh_p1 / enh_p2)
        other_Enr_list.append(other_p1 / other_p2)

        # Calculate the standard errors of log odds ratios and enrichment values.
        SE_logOR_pro = sqrt((1 / float(seg_method_dict[seg_method][0][0])) + (
            1 / (float(effective_mut_num) - float(seg_method_dict[seg_method][0][0]))) + (
                                1 / float(seg_method_dict[seg_method][1][0])) + (
                                1 / (float(chromosome_res_num) - float(seg_method_dict[seg_method][1][0]))))
        SE_logOR_enh = sqrt((1 / float(seg_method_dict[seg_method][0][1])) + (
            1 / (float(effective_mut_num) - float(seg_method_dict[seg_method][0][1]))) + (
                                1 / float(seg_method_dict[seg_method][1][1])) + (
                                1 / (float(chromosome_res_num) - float(seg_method_dict[seg_method][1][1]))))
        SE_logOR_other = sqrt((1 / float(seg_method_dict[seg_method][0][2])) + (
            1 / (float(effective_mut_num) - float(seg_method_dict[seg_method][0][2]))) + (
                                  1 / float(seg_method_dict[seg_method][1][2])) + (
                                  1 / (float(chromosome_res_num) - float(seg_method_dict[seg_method][1][2]))))
        SE_logEnr_pro = sqrt(float((effective_mut_num - seg_method_dict[seg_method][0][0])) / (
            effective_mut_num * seg_method_dict[seg_method][0][0]) + float(
            (chromosome_res_num - seg_method_dict[seg_method][1][0])) / (
                                 chromosome_res_num * seg_method_dict[seg_method][1][0]))
        SE_logEnr_enh = sqrt(float((effective_mut_num - seg_method_dict[seg_method][0][1])) / (
            effective_mut_num * seg_method_dict[seg_method][0][1]) + float(
            (chromosome_res_num - seg_method_dict[seg_method][1][1])) / (
                                 chromosome_res_num * seg_method_dict[seg_method][1][1]))
        SE_logEnr_other = sqrt(float((effective_mut_num - seg_method_dict[seg_method][0][2])) / (
            effective_mut_num * seg_method_dict[seg_method][0][2]) + float(
            (chromosome_res_num - seg_method_dict[seg_method][1][2])) / (
                                   chromosome_res_num * seg_method_dict[seg_method][1][2]))
        SE_OR_pro_list.append(SE_logOR_pro * (pro_p1 / (1 - pro_p1)) / (pro_p2 / (1 - pro_p2)))
        SE_OR_enh_list.append(SE_logOR_enh * (enh_p1 / (1 - enh_p1)) / (enh_p2 / (1 - enh_p2)))
        SE_OR_other_list.append(SE_logOR_other * (other_p1 / (1 - other_p1)) / (other_p2 / (1 - other_p2)))
        SE_Enr_pro_list.append(SE_logEnr_pro * pro_p1 / pro_p2)
        SE_Enr_enh_list.append(SE_logEnr_enh * enh_p1 / enh_p2)
        SE_Enr_other_list.append(SE_logEnr_other * other_p1 / other_p2)

        Z_pro = log(pro_p1 / pro_p2) / SE_logEnr_pro
        Z_enh = log(enh_p1 / enh_p2) / SE_logEnr_enh
        Z_other = log(other_p1 / other_p2) / SE_logEnr_other
        Z_list.append([Z_pro, Z_enh, Z_other])

    # Use matplotlib to plot the enrichment with error bars of the mutations in promoters, enhancers and other regions.
    N = len(seg_method_dict)
    ind = np.arange(N)
    width = 0.25
    fig, ax = plt.subplots()
    pro_plot = [x - 1 for x in pro_Enr_list]
    enh_plot = [x - 1 for x in enh_Enr_list]
    other_plot = [x - 1 for x in other_Enr_list]
    rects1 = ax.bar(ind + 0.1, pro_plot, width, color='#fbb4ae', yerr=SE_Enr_pro_list)
    rects2 = ax.bar(ind + width + 0.1, enh_plot, width, color='#b3cde3', yerr=SE_Enr_enh_list)
    rects3 = ax.bar(ind + 2 * width + 0.1, other_plot, width, color='#ccebc5', yerr=SE_Enr_other_list)

    ax.axhline(0, color='black', lw=1)
    plt.gca().yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + 1))

    # Mark the statistically significant ones with astericks.
    for i in range(len(seg_method_list)):
        for j in range(3):
            Z_score = Z_list[i][j]
            if abs(Z_score) > 1.96:
                if abs(Z_score) > 2.577:
                    if abs(Z_score) > 3.3:
                        mark = '***'
                    else:
                        mark = '**'
                else:
                    mark = '*'
            else:
                mark = ''
            if j == 0:  # Promoter
                ax.annotate(mark, xy=(0.1 + ind[i] + 0.5 * width, pro_plot[i] + SE_Enr_pro_list[i] + 0.03),
                            zorder=10, ha='center')
            elif j == 1:  # Enhancer
                ax.annotate(mark, xy=(0.1 + ind[i] + 1.5 * width, enh_plot[i] + SE_Enr_enh_list[i] + 0.03),
                            zorder=10, ha='center')
            else:  # Other
                ax.annotate(mark, xy=(0.1 + ind[i] + 2.5 * width, other_plot[i] + SE_Enr_other_list[i] + 0.03),
                            zorder=10, ha='center')
    
    ax.set_ylim(-0.5, 1.0)
    ax.legend((rects1[0], rects2[0], rects3[0]), ('TSS', 'Enhancer', 'Other'), loc=2)
    ax.set_ylabel('Enrichment of HGMD non-coding mutations')
    ax.set_xticks(ind + 1.5 * width + 0.1)
    ax.set_xticklabels(seg_method_list)
    plt.savefig("overall_only_strong.png", dpi=1000)

# Now we combine the DM, DM?, DP and DFP mutations.
valid_mut_dict = {}
for iter in range(len(DM_dict)):
    valid_mut_dict[iter] = DM_dict[iter]
for iter in range(len(DMq_dict)):
    valid_mut_dict[iter + len(DM_dict)] = DMq_dict[iter]
for iter in range(len(DP_dict)):
    valid_mut_dict[iter + len(DM_dict) + len(DMq_dict)] = DP_dict[iter]
for iter in range(len(DFP_dict)):
    valid_mut_dict[iter + len(DM_dict) + len(DMq_dict) + len(DP_dict)] = DFP_dict[iter]

calculate_overall_enrichment(valid_mut_dict)
