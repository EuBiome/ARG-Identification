import operator
import matplotlib.pyplot as plt
import matplotlib
import argparse
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

from constants import SPECIES, GENRE

matplotlib.style.use('ggplot')
from pandas import DataFrame
from helpers import find_taxa_path, cnt_per_drug_abundance, extract_species_name, link_taxa_antibiotic, group_by_taxa, \
    to_proper_form, prepare_input_for_graph, draw_graph

# This is a sample Python script.

# Press Maiusc+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

#
# def print_hi(name):
#     # Use a breakpoint in the code line below to debug your script.
#     print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
#
#
# # Press the green button in the gutter to run the script.
# if __name__ == '__main__':
#     print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/


# import subprocess
#
# subprocess.call(['sh', '/home/margherita/Scrivania/Eubiome/pipe_mod.bash'])

ap = argparse.ArgumentParser()

ap.add_argument("--level_of_interest", required=True, help="Write S if you want to produce the graph at species level and G if you want it at genre level")
ap.add_argument("--rgi_output", required=True, help="Insert the path to the RGI output")
ap.add_argument("--kraken2_output", required=True, help="path to the kraken2 output")
ap.add_argument("--kraken2_report", required=True, help="path to the kraken2 report")
#ap.add_argument("--pathogens-list", required=True, help="path of the file containg the pathogens list")
ap.add_argument("--rgi_analysis", required=True,
                help="Insert the path where you want to store the result of your analysis")

args = vars(ap.parse_args())
level_of_interest = args["level_of_interest"]
rgi_output = args["rgi_output"]
kraken2_output = args["kraken2_output"]
rgi_analysis = args["rgi_analysis"]
#pathogens = args["pathogens_list"]
kraken2_report = args["kraken2_report"]

with open(rgi_output, mode='r') as rgi_file, open(kraken2_output, mode='r') as k2_results_file, open(kraken2_report, mode="r") as report_file:
    # list of all the contigs contigs associated to the ORFs identified by rgi
    orf_contigs = []
    # [['contig_id1', ['macrolide antibiotic', 'fluoroquinolone antibiotic',...]], ['contig_id2',['macrolide
    # antibiotic', 'lincosamide antibiotic',..]],...]
    info_array = []
    # List of all the antibiotic classes for by rgi
    drug_lists = []
    cnt = 0

    for line in rgi_file:
        cnt += 1
        if cnt >= 2:
            info_element = []
            id_contig = line.split()[0]
            id_contig = id_contig[:id_contig.rfind("_")]
            orf_contigs.append(id_contig)
            contig_drugs = line.split("\t")[14]
            contig_drugs_list = []
            for j in range(len(contig_drugs.split("; "))):
                drug_class = contig_drugs.split("; ")[j]
                contig_drugs_list.append(drug_class)
            info_element.append(id_contig)
            info_element.append(contig_drugs_list)
            drug_lists.extend(contig_drugs_list)
            info_array.append(info_element)

    # Delete duplicates and order by name
    drug_lists = sorted(drug_lists, key=operator.itemgetter(0), reverse=True)
    drug_lists = list(set(drug_lists))

    orf_contigs = list(set(orf_contigs))

    # contigs of interest associated with their labels
    labelled_contigs = []
    # filtering of the k2 by keeping only contigs associated with orfs identified by rgi
    for line2 in k2_results_file:
        if line2.split()[1] in orf_contigs:
            labbeled_contig = [line2.split()[1], line2.split("\t")[2].split(" (")[0],
                               line2.split("taxid ")[1].split(")")[0]]
            labelled_contigs.append(labbeled_contig)
            # print(labbeled_contig)

    # Trasformo il file del report in una lista di liste:
    # Ogni sottolista è così strutturata: ex. ['G', '10494', 'Lymphocystivirus']
    # Perciò: (livello tassonomico, id_tassonomico, nomenclatura_tassonomica)

    report_to_lists = []
    for line3 in report_file:
        report_to_list = [line3.split()[3], line3.split()[4], line3.split("\t")[5].lstrip()[:-1]]
        report_to_lists.append(report_to_list)

    def classification_level(specie, genre, family):
        if specie != "---":
            level = "S"
        else:
            if genre != "---":
                level = "G"
            else:
                if family != "---":
                    level = "F"
                else:
                    level = "U"
        return level

    find_taxa_path(labelled_contigs, report_to_lists)

    for i in range(len(labelled_contigs)):
        labelled_contigs[i][3] = extract_species_name(labelled_contigs[i][3], labelled_contigs[i][4])

    drug_abundance_list = cnt_per_drug_abundance(info_array, drug_lists)
    drug_abundance_list = sorted(drug_abundance_list, key=operator.itemgetter(0), reverse=True)
    drug_abundance_list = [tuple(l) for l in drug_abundance_list]
    drug_abundance_list = dict(drug_abundance_list)

    all_drug_keys = frozenset(drug_abundance_list)

    # [('triclosan', 1), ('tetracycline antibiotic', 42), ('streptogramin antibiotic', 19), ... ]
    drug_abundance_list = list(drug_abundance_list.items())

    # list of lists: [[tax_id, [list of antibiotic]], ...]
    if level_of_interest.casefold() == "s":
        taxa_antibiotic_associators = link_taxa_antibiotic(info_array, labelled_contigs, SPECIES)
    else:
        taxa_antibiotic_associators = link_taxa_antibiotic(info_array, labelled_contigs, GENRE)

    #taxa_antibiotic_associators_genre = link_taxa_antibiotic(info_array, labelled_contigs, GENRE)

    # sorting the list according to the taxid
    sorted_taxa_antibiotic_associators = sorted(taxa_antibiotic_associators, key=operator.itemgetter(0))
    #sorted_taxa_antibiotic_associators_genre = sorted(taxa_antibiotic_associators_genre, key=operator.itemgetter(0))

    # if two tax-id are the same I merge the antibiotic associated to that tax-id
    # merged_taxa_antibiotic_associators_genre = group_by_taxa(sorted_taxa_antibiotic_associators_genre)
    merged_taxa_antibiotic_associators = group_by_taxa(sorted_taxa_antibiotic_associators)

    #for i in merged_taxa_antibiotic_associator_species: print(i)

    final_taxa_antibiotic_associators = to_proper_form(merged_taxa_antibiotic_associators, all_drug_keys,
                                                       drug_abundance_list)

    # [['Clostridioides Difficile', 0, 0, 0, 0], .. ['Clostridium Perfringens', 1, 0, 1, 0]]
    # Per each taxa, quantity of antibiotic resistance found per each antibiotic
    taxa_quantity_associators = prepare_input_for_graph(final_taxa_antibiotic_associators)[1]
    # antibiotic list = list of all the antibiotic indentified without repetitions
    antibiotic_list = prepare_input_for_graph(final_taxa_antibiotic_associators)[0]

    rgi_analysis_level = rgi_analysis + "_" + level_of_interest
    draw_graph(taxa_quantity_associators, sns, antibiotic_list, rgi_analysis_level)
