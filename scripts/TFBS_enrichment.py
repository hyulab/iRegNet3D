# Codes used for the analyses in Fig 4a.

from math import sqrt
from math import log
import interval

def measure(x):
    from interval import fpu
    return fpu.up(lambda: sum((c.sup - c.inf for c in x), 0))

# Read in the TF motif file with a specified threshold.
mut_TF_reg_dict = {}
def read_mut_file_with_threshold(threshold = 6):
    mut_TF_reg_dict.clear()
    TF_interval_dict = {}
    search_interval_dict = {}
    with open("mut_scan.bed.txt") as mut_file:
        mut_reader = mut_file.readlines()
        for row in mut_reader:
            entry = row.split()
            if float(entry[4]) > threshold:
                chr = entry[0][3:]
                start = int(entry[1])
                end = int(entry[2])
                TF_name = entry[3]
                TF_interval = interval.interval[start, end]
                search_interval = interval.interval[int(entry[6]), int(entry[7])]
                if chr not in mut_TF_reg_dict:
                    mut_TF_reg_dict[chr] = [(start, end, TF_name)]
                    TF_interval_dict[chr] = [TF_interval]
                    search_interval_dict[chr] = [search_interval]
                else:
                    mut_TF_reg_dict[chr].append((start, end, TF_name))
                    TF_interval_dict[chr].append(TF_interval)
                    search_interval_dict[chr].append(search_interval)
    mut_file.close()
    # Calculate p2.
    TF_residues = 0
    search_residues = 0
    for chr in TF_interval_dict:
        TF_intervals_in_chr = interval.interval.union(TF_interval_dict[chr])
        search_intervals_in_chr = interval.interval.union(search_interval_dict[chr])
        TF_residues += measure(TF_intervals_in_chr)
        search_residues += measure(search_intervals_in_chr)
    return [float(TF_residues) / search_residues, TF_residues, search_residues]

#############################################################################################################
# For a series of thresholds, calculate the percentage of non-coding mutations inside TF binding motifs.

thresholds = [6, 7, 8, 9]
mut_results = []
mut_enrichment = []
mut_SE_enrichment = []
for threshold in thresholds:
    in_count = 0
    returned = read_mut_file_with_threshold(threshold)
    p2 = returned[0]
    TF_res = returned[1]
    all_res = returned[2]
    for mut_num in valid_mut_dict:
        mut_chr = valid_mut_dict[mut_num][0]
        mut_pos = valid_mut_dict[mut_num][1]
        for TF_region in mut_TF_reg_dict[mut_chr]:
            if TF_region[0] <= mut_pos <= TF_region[1]:
                in_count += 1
                break
    p1 = float(in_count)/1666
    mut_results.append(p1)
    SE_enrichment = sqrt(1/float(in_count) - 1/float(1666) + 1/float(TF_res) - 1/float(all_res)) * p1 / p2
    mut_enrichment.append(p1 / p2)
    mut_SE_enrichment.append(SE_enrichment)
    print in_count, TF_res, all_res
thresholds = [6, 7, 8, 9]

#############################################################################################################
# For the three population SNP groups, calculate the percentage of them inside TF binding motifs.
def read_SNP_file(file_name):
    SNP_dict = {}
    with open(file_name) as SNP_file:
        for row in SNP_file:
            entry = row.split()
            chr = entry[0][3:]
            pos = int(entry[1])
            if chr not in SNP_dict:
                SNP_dict[chr] = [pos]
            else:
                SNP_dict[chr].append(pos)
    return SNP_dict

less_than_1_dict = read_SNP_file("less_than_1.bed")
one_to_10_dict = read_SNP_file("one_to_10.bed")
over_10_dict = read_SNP_file("over_10.bed")

def read_scan_result(scan_file_name, threshold=6):
    SNP_TF_reg_dict = {}
    TF_interval_dict = {}
    search_interval_dict = {}
    with open(scan_file_name) as scan_file:
        scan_reader = scan_file.readlines()
        for row in scan_reader:
            scan_entry = row.split()
            if float(scan_entry[4]) > threshold:
                chr = scan_entry[0][3:]
                start = int(scan_entry[1])
                end = int(scan_entry[2])
                TF_name = scan_entry[3]
                TF_interval = interval.interval[start, end]
                search_interval = interval.interval[int(scan_entry[6]), int(scan_entry[7])]
                if chr not in SNP_TF_reg_dict:
                    SNP_TF_reg_dict[chr] = [(start, end, TF_name)]
                    TF_interval_dict[chr] = [TF_interval]
                    search_interval_dict[chr] = [search_interval]
                else:
                    SNP_TF_reg_dict[chr].append((start, end, TF_name))
                    TF_interval_dict[chr].append(TF_interval)
                    search_interval_dict[chr].append(search_interval)
    # Calculate p2.
    TF_residues = 0
    search_residues = 0
    for chr in TF_interval_dict:
        TF_intervals_in_chr = interval.interval.union(TF_interval_dict[chr])
        search_intervals_in_chr = interval.interval.union(search_interval_dict[chr])
        TF_residues += measure(TF_intervals_in_chr)
        search_residues += measure(search_intervals_in_chr)
    return [SNP_TF_reg_dict, TF_residues, search_residues]

def calculate_result(SNP_dict, scan_file_name):
    result_list = []
    enrichment_list = []
    SE_Enr_list = []
    for threshold in thresholds:
        scan_returned = read_scan_result(scan_file_name, threshold)
        scan_dict = scan_returned[0]
        TF_res = scan_returned[1]
        all_res = scan_returned[2]
        in_count = 0
        for chr in SNP_dict:
            for pos in SNP_dict[chr]:
                for TF_region in scan_dict[chr]:
                    if TF_region[0] <= pos <= TF_region[1]:
                        in_count += 1
                        break
        p1 = float(in_count)/2000
        result_list.append(p1)
        p2 = float(TF_res) / all_res
        enrichment_list.append(p1/p2)
        SE_Enr = sqrt(1/float(in_count) - 1/float(2000) + 1/float(TF_res) - 1/float(all_res)) * p1 / p2
        SE_Enr_list.append(SE_Enr)
        print in_count, TF_res, all_res
    return result_list, enrichment_list, SE_Enr_list

less_than_1_results, less_than_1_enr, less_than_1_SE_enr = calculate_result(less_than_1_dict, "less_than_1_scan.bed")
one_to_10_results, one_to_10_enr, one_to_10_SE_enr = calculate_result(one_to_10_dict, "one_to_10_scan.bed")
over_10_results, over_10_enr, over_10_SE_enr = calculate_result(over_10_dict, "over_10_scan.bed")

