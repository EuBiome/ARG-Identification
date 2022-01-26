import argparse
from helpers import k2_output_to_matrix, compute_score, is_C, arg_to_matrix

ap = argparse.ArgumentParser()

ap.add_argument("--plasmid_classification_path", required=True, help="Insert the path of plasmid classification")
ap.add_argument("--microbial_classification_path", required=True, help="Insert the path of the microbial classification")
ap.add_argument("--microbial_report_path", required=True, help="Insert the path of the microbial report")
ap.add_argument("--rgi_path", required=True, help="Insert the path of the rgi output")
ap.add_argument("--genes_list", required=True, help="Insert the path of the genes list")
ap.add_argument("--pathogens_list", required=True, help="Insert the path of the pathogens list")
args = vars(ap.parse_args())

plasmid_classification_path = args["plasmid_classification_path"]
microbial_classification_path = args["microbial_classification_path"]
microbial_report_path = args["microbial_report_path"]
rgi_path = args["rgi_path"]
genes_list = args["genes_list"]
pathogens_list = args["pathogens_list"]


path_len = len(microbial_classification_path.split("/"))
base_path_pieces = microbial_classification_path.split("/")[0:path_len-1]
base_path = ""

for i in range(len(base_path_pieces)):
    base_path = base_path + base_path_pieces[i] + "/"

lista_1 = base_path + "lista1"
lista_2 = base_path + "lista2"
lista_3 = base_path + "lista3"
lista_4 = base_path + "lista4"

with open(rgi_path, mode="r") as rgi_file, open(plasmid_classification_path, mode="r") as plasmid_classification_file,\
        open(microbial_classification_path, mode="r") as microbial_classification_file, \
        open(microbial_report_path, mode="r") as microbial_report_file, open(genes_list, mode="r") as genes_list, open(lista_1, mode="w") as list_1, open(lista_2, mode="w") as list_2, \
        open(lista_3, mode="w") as list_3, open(lista_4, mode="w") as list_4, open(pathogens_list, mode="r") as pathogens_file:

    # Trasformo il file di rgi in una lista di liste. Il primo elemento della lista è la contig-id e quelli successivi
    # nomi dei arg presenti al suo interno es: ['k127_82339', "AAC(6')-Im;APH(2'')-IIa"]
    rgi_matrix = []
    cnt1 = 0
    for rgi_line in rgi_file:
        if cnt1 > 0:
            contig_des = rgi_line.split("\t")[1]
            contig = str(contig_des.split("_")[0] + "_" + contig_des.split("_")[1])
            arg = str(rgi_line.split("\t")[8])
            if cnt1 == 1:
                rgi_array = []
                rgi_array.append(contig)
                rgi_array.append(arg)
                rgi_matrix.append(rgi_array)
            elif cnt1 > 1:
                if contig == rgi_matrix[len(rgi_matrix) - 1][0]:
                    rgi_matrix[len(rgi_matrix) - 1][1] = str(rgi_matrix[len(rgi_matrix) - 1][1] + ';' + arg)
                else:
                    rgi_array = []
                    rgi_array.append(contig)
                    rgi_array.append(arg)
                    rgi_matrix.append(rgi_array)
        cnt1 += 1

    empty_plasmid_matrix = []
    empty_minikraken_matrix = []

    # trasformo i file di output di kraken2 in liste di liste. es: ogni lista è così composta: ['U', 'k127_138004', 'unclassified (taxid 0)', '0:327']
    plasmid_matrix = k2_output_to_matrix(plasmid_classification_file, empty_plasmid_matrix)
    minikraken_matrix = k2_output_to_matrix(microbial_classification_file, empty_minikraken_matrix)

    both_classified = 0
    plasmid_classified = 0
    minikraken_classified = 0
    unclassified = 0
    # Per ogni contig c'è scritto se è stata classificata, la classificazione e gli score per i rispettivi database
    # list of list. ex list: ['U', 'k127_138049', 'unclassified (taxid 0)', 1.0, 'unclassified (taxid 0)', 1.0]
    classification_matrix = []

    for contig in range(len(plasmid_matrix)):
        array = []
        is_C_plasmid = plasmid_matrix[contig][0]
        is_C_microbial = minikraken_matrix[contig][0]
        is_C_string = is_C(is_C_plasmid, is_C_microbial)
        array.append(is_C_string)
        array.append(plasmid_matrix[contig][1])
        array.append(plasmid_matrix[contig][2])
        plasmid_assigned_taxid = plasmid_matrix[contig][2].split("taxid ")[1][:-1]
        score_plasmid = compute_score(plasmid_matrix[contig][3], plasmid_assigned_taxid)
        array.append(score_plasmid)
        array.append(minikraken_matrix[contig][2])
        minikraken_assigned_taxid = minikraken_matrix[contig][2].split("taxid ")[1][:-1]
        score_minikraken = compute_score(minikraken_matrix[contig][3], minikraken_assigned_taxid)
        array.append(score_minikraken)
        classification_matrix.append(array)
        if plasmid_matrix[contig][0] == "C" and minikraken_matrix[contig][0] == "C":
            both_classified += 1
        elif plasmid_matrix[contig][0] == "C" and minikraken_matrix[contig][0] == "U":
            plasmid_classified += 1
        elif plasmid_matrix[contig][0] == "U" and minikraken_matrix[contig][0] == "C":
            minikraken_classified += 1
        else:
            unclassified += 1

    # Aggiungo le informazioni relative ai geni di antibiotico resistenza rilevati
    for classification_list in classification_matrix:
        contig = classification_list[1]
        for arg_list in rgi_matrix:
            if arg_list[0] == contig:
                classification_list.append(arg_list[1])

    arg_matrix = arg_to_matrix(genes_list)
    #print(arg_matrix)

    # Genero il file 1
    for line in classification_matrix:
        if len(line)>6:
            #print(line)
            identified_genes = line[6]
            number_of_genes = len(identified_genes.split(";"))
            if number_of_genes == 1:
                for list in arg_matrix:
                    genes_names_number = len(list)
                    for i in range(genes_names_number):
                        gene_name = list[i]
                        if identified_genes == gene_name:
                            list_1.write(list[0] + "\t" + line[4] + "\t" + line[1] + "\n")
                            line.append("bad_gene")
            else:
                for j in range(number_of_genes):
                    identified_gene = identified_genes.split(";")[j]
                    #print(identified_gene)
                    for list in arg_matrix:
                        genes_names_number = len(list)
                        for i in range(genes_names_number):
                            gene_name = list[i]
                            if identified_gene == gene_name:
                                list_1.write(list[0] + "\t" + line[4] + "\t" + line[1] + "\n")
                                if len(line)< 8: line.append("bad_gene")

    
    pathogens_set = set()
    for line in pathogens_file:
        pathogens_set.add(line[:-1])
    
    for line in classification_matrix:
        if len(line) == 7:
            if line[2] == "unclassified (taxid 0)":
                list_4.write(line[6] + "\t" + line[4] + "\t" + line[1] + "\n")
            else:
                if line[4].split("(taxid ")[1].split(")")[0] in pathogens_set:
                    list_2.write(line[6] + "\t" + line[4] + "\t" + line[1] + "\n")
                else:
                    list_3.write(line[6] + "\t" + line[4] + "\t" + line[1] + "\n")
