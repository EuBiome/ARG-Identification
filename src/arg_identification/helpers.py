def compute_score(kmer_groups, assigned_taxid):
    C = 0
    Q = 0
    for i in range(len(kmer_groups.split())):
        taxid = kmer_groups.split()[i].split(":")[0]
        kmer_num = kmer_groups.split()[i].split(":")[1]
        Q = Q + float(kmer_num)
        if taxid == assigned_taxid:
            C = C + int(kmer_num)
    score = C/Q
    return score

def is_C(is_C_plasmid, is_C_microbial):
    if is_C_microbial == "C" or is_C_plasmid == "C":
        is_C_string = "C"
    else:
        is_C_string = "U"
    return is_C_string

def k2_output_to_matrix(file, empty_matrix):
    for line in file:
        line_array = []
        line_array.append(line.split("\t")[0])
        line_array.append(line.split("\t")[1])
        line_array.append(line.split("\t")[2])
        line_array.append(line.split("\t")[4][:-1])
        empty_matrix.append(line_array)
    return empty_matrix

def arg_to_matrix(genes_list):
    genes_matrix = []
    for line in genes_list:
        genes_array = line[:-1].split(";")
        genes_matrix.append(genes_array)
    return genes_matrix
