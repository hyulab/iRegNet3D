# Helper script for randomly generating a sample of population SNPs. 

less_than_1_dict = {}
one_to_10_dict = {}
over_10_dict = {}

less_than_1_iter = 0
one_to_10_iter = 0
over_10_iter = 0

with open("snp146Common.txt") as snp_file:
    for row in snp_file:
        entry = row.split()
        genomic = entry[10]
        category = entry[11]
        allele_f = entry[23]
        allele_fs = allele_f.split(',')
        MAF = min(float(allele_fs[0]), float(allele_fs[1]))
        chrom = entry[1][3:]
        pos = int(entry[2])
        if genomic == 'genomic' and category == 'single' and chrom in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X']:
            if MAF < 0.01:
                less_than_1_dict[less_than_1_iter] = [chrom, pos, MAF]
                less_than_1_iter += 1
            elif MAF >= 0.01 and MAF < 0.1:
                one_to_10_dict[one_to_10_iter] = [chrom, pos, MAF]
                one_to_10_iter += 1
            elif MAF >= 0.1 and MAF < 0.5:
                over_10_dict[over_10_iter] = [chrom, pos, MAF]
                over_10_iter += 1
    snp_file.close()

import random
foo = random.SystemRandom()
set_count = 2000

def collect(factor, limit, my_dict, output_file):
    with open(output_file, 'w+') as out_file:
        counter = 0
        while(counter < set_count):
            rand_num = int(factor * foo.random())
            if rand_num < limit:
                counter += 1
                entry = my_dict[rand_num]
                new_line = 'chr' + entry[0] + '\t' + str(entry[1] - 2000) + '\t' + str(entry[1] + 2000) + '\n'
                out_file.write(new_line)

collect(100000, len(less_than_1_dict), less_than_1_dict, 'less_than_1.bed')
collect(100000, len(one_to_10_dict), one_to_10_dict, 'one_to_10.bed')
collect(100000, len(over_10_dict), over_10_dict, 'over_10.bed')
            
