import operator
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.style.use('ggplot')
from pandas import DataFrame


def extract_species_name(extended_name, genre):
    if extended_name != "---" and genre != "---":
        short_name = extended_name.replace(genre, "")
        if " sp. " in short_name:
            short_name = short_name.replace(" sp. ", "")
        else:
            short_name = short_name[1:]
        if " subsp." in short_name:
            short_name = short_name.split(" subsp.")[0]
    elif extended_name != "---" and genre == "---":
        short_name = extended_name
        if "sp. " in short_name:
            short_name = short_name.replace("sp. ", "")
        if " subsp." in short_name:
            short_name = short_name.split(" subsp.")[0]
    else:
        short_name = extended_name
    return short_name


# For each antibiotic I count how many antibiotic resistance gene are associated to it
def cnt_per_drug_abundance(info_array, drug_lists):
    drug_abundance_list = []
    for k in range(len(drug_lists)):
        drug_abundance = [drug_lists[k], int(0)]
        drug_abundance_list.append(drug_abundance)

    for i in range(len(info_array)):
        for j in range(len(info_array[i][1])):
            j_index = drug_lists.index(info_array[i][1][j])
            drug_abundance_list[j_index][1] = drug_abundance_list[j_index][1] + 1
    return drug_abundance_list


def link_taxa_antibiotic(info_array, labelled_contigs, level):
    taxa_antibiotic_associators = []
    for i in range(len(info_array)):
        for j in range(len(labelled_contigs)):
            # Se stiamo trattando la stessa contig
            if labelled_contigs[j][0] == info_array[i][0]:
                if level == "species":
                    if labelled_contigs[j][3] != "---":
                        taxa_antibiotic_associator = [labelled_contigs[j][4] + " " + labelled_contigs[j][3],
                                                      info_array[i][1]]
                        taxa_antibiotic_associators.append(taxa_antibiotic_associator)
                        break
                    else:
                        if labelled_contigs[j][4] != "---":
                            taxa_antibiotic_associator = [labelled_contigs[j][4] + " (" + labelled_contigs[j][5] + " )", info_array[i][1]]
                            taxa_antibiotic_associators.append(taxa_antibiotic_associator)
                        else:
                            taxa_antibiotic_associator = [labelled_contigs[j][4] + " " + labelled_contigs[j][3],
                                                          info_array[i][1]]
                            taxa_antibiotic_associators.append(taxa_antibiotic_associator)
                else:
                    taxa_antibiotic_associator = [labelled_contigs[j][5] + " " + labelled_contigs[j][4],
                                                  info_array[i][1]]
                    taxa_antibiotic_associators.append(taxa_antibiotic_associator)
    return taxa_antibiotic_associators


def group_by_taxa(sorted_taxa_antibiotic_associators):
    prev_tax_id = "N"
    merged_taxa_antibiotic_associators = []
    for i in range(len(sorted_taxa_antibiotic_associators)):
        current_tax_id = sorted_taxa_antibiotic_associators[i][0]
        if prev_tax_id != current_tax_id:
            new_list = [current_tax_id, sorted_taxa_antibiotic_associators[i][1]]
            merged_taxa_antibiotic_associators.append(new_list)
            prev_tax_id = current_tax_id
        else:
            merged_taxa_antibiotic_associators[(len(merged_taxa_antibiotic_associators)) - 1][1].extend(
                sorted_taxa_antibiotic_associators[i][1])
    return merged_taxa_antibiotic_associators


def to_proper_form(merged_taxa_antibiotic_associators, all_drug_keys, drug_abundance_list):
    final_taxa_antibiotic_associators = []
    # Creo una lista di liste dove il primo elemento è un tax-id ed il secondo un dizionario. Ciascun dizionario ha come
    # chiave l'antibiotico e come valore il numero di resistenze associate a quell'antibiotico trovate per quel tax-id.
    for j in range(len(merged_taxa_antibiotic_associators)):
        final_taxa_antibiotic_associator = []
        my_dict = {i: merged_taxa_antibiotic_associators[j][1].count(i) for i in
                   merged_taxa_antibiotic_associators[j][1]}
        final_taxa_antibiotic_associator.append(merged_taxa_antibiotic_associators[j][0])
        final_taxa_antibiotic_associator.append(my_dict)
        final_taxa_antibiotic_associators.append(final_taxa_antibiotic_associator)

    # Aggiungo ai dizionari le chiavi con valori pari a zero
    for i in range(len(final_taxa_antibiotic_associators)):
        for drug in all_drug_keys:
            if drug not in final_taxa_antibiotic_associators[i][1]:
                final_taxa_antibiotic_associators[i][1][drug] = 0

    # Trasformo i dizionari in liste di liste e ordino le sottoliste in ordine alfabetico rispetto al nome del antibiotico
    for i in range(len(final_taxa_antibiotic_associators)):
        final_taxa_antibiotic_associators[i][1] = sorted(list(final_taxa_antibiotic_associators[i][1].items()),
                                                         key=operator.itemgetter(0), reverse=True)
        final_taxa_antibiotic_associators[i][1] = [list(ele) for ele in final_taxa_antibiotic_associators[i][1]]

    # Riordino questa volta rispetto al numero totale di antibiotico resistenze trovate per una dato antibiotico
    for i in range(len(final_taxa_antibiotic_associators)):
        for j in range(len(final_taxa_antibiotic_associators[i][1])):
            final_taxa_antibiotic_associators[i][1][j].append(drug_abundance_list[j][1])
        final_taxa_antibiotic_associators[i][1] = sorted(final_taxa_antibiotic_associators[i][1],
                                                         key=operator.itemgetter(2), reverse=True)
    return final_taxa_antibiotic_associators


def prepare_input_for_graph(final_taxa_antibiotic_associators):
    antibiotic_list = []
    taxa_quantity_associators = []
    for i in range(len(final_taxa_antibiotic_associators)):
        taxa_quantity_associator = [final_taxa_antibiotic_associators[i][0]]
        if i == 0:
            for j in range(len(final_taxa_antibiotic_associators[i][1])):
                antibiotic_list.append(final_taxa_antibiotic_associators[i][1][j][0])
        for j in range(len(final_taxa_antibiotic_associators[i][1])):
            taxa_quantity_associator.append(final_taxa_antibiotic_associators[i][1][j][1])
        taxa_quantity_associators.append(taxa_quantity_associator)
    return antibiotic_list, taxa_quantity_associators


def draw_graph(taxa_quantity_associators, sns, antibiotic_list, rgi_analysis):
    data_dict = {item[0]: item[1:] for item in taxa_quantity_associators}
    sns.set_theme()
    colors = sns.color_palette("cubehelix", n_colors=len(taxa_quantity_associators))
    cmap1 = LinearSegmentedColormap.from_list("my_colormap", colors)
    index = antibiotic_list
    df = DataFrame(data_dict, index=index)
    ax = df.plot(kind='bar', stacked=True, colormap=cmap1, figsize=(11, 8))
    ax.legend(loc=(1.01, 0.01), ncol=2)
    plt.tight_layout()
    plt.savefig(rgi_analysis)


def find_taxa_path(labelled_contigs, rep_to_lists):
    for i in range(len(labelled_contigs)):
        taxid_any_level = labelled_contigs[i][2]
        for j in range(len(rep_to_lists)):
            if taxid_any_level == rep_to_lists[j][1]:
                actual_tax_level = rep_to_lists[j][0]
                ####################################################################################################
                if actual_tax_level.startswith("S") and len(actual_tax_level) > 1:
                    # sono a livello sottospecie
                    n = 1
                    while rep_to_lists[j - n][0] != "S":
                        n += 1
                    labelled_contigs[i].append(rep_to_lists[j - n][2])
                    quadra = False
                    if "[" in rep_to_lists[j - n][2]:
                        quadra = True
                    n = 1
                    while rep_to_lists[j - n][0] != "G":
                        if rep_to_lists[j - n][0].startswith("S") == False and rep_to_lists[j - n][0].startswith(
                                "G") == False:
                            labelled_contigs[i][3] = "---"
                            labelled_contigs[i].append("---")
                            break
                        n += 1
                    if len(labelled_contigs[i]) < 5:
                        if not quadra:
                            labelled_contigs[i].append(rep_to_lists[j - n][2])
                        else:
                            labelled_contigs[i].append("[" + rep_to_lists[j - n][2] + "]")

                    n = 1
                    while rep_to_lists[j - n][0] != "F":
                        if rep_to_lists[j - n][0].startswith("S") == False and rep_to_lists[j - n][0].startswith(
                                "G") == False and rep_to_lists[j - n][0].startswith("F") == False:
                            labelled_contigs[i].append("Undefined")
                            break
                        n += 1
                    if len(labelled_contigs[i]) < 6:
                        labelled_contigs[i].append(rep_to_lists[j - n][2])
                ####################################################################################################
                elif actual_tax_level == "S":
                    labelled_contigs[i].append(labelled_contigs[i][1])
                    quadra = False
                    if "[" in labelled_contigs[i][1]:
                        quadra = True
                    n = 1
                    while rep_to_lists[j - n][0] != "G":
                        if rep_to_lists[j - n][0].startswith("S") == False and rep_to_lists[j - n][0].startswith(
                                "G") == False:
                            labelled_contigs[i][3] = "---"
                            labelled_contigs[i].append("---")
                            break
                        n += 1

                    if len(labelled_contigs[i]) < 5:
                        if not quadra:
                            labelled_contigs[i].append(rep_to_lists[j - n][2])
                        else:
                            labelled_contigs[i].append("[" + rep_to_lists[j - n][2] + "]")

                    n = 1
                    while rep_to_lists[j - n][0] != "F":
                        if rep_to_lists[j - n][0].startswith("S") == False and rep_to_lists[j - n][0].startswith(
                                "G") == False and rep_to_lists[j - n][0].startswith("F") == False:
                            labelled_contigs[i].append("Undefined")
                            break
                        n += 1
                    if len(labelled_contigs[i]) < 6:
                        labelled_contigs[i].append(rep_to_lists[j - n][2])
                ####################################################################################################
                elif actual_tax_level.startswith("G") and len(actual_tax_level) > 1:
                    # sono a livello sottospecie
                    labelled_contigs[i].append("---")
                    n = 1
                    while rep_to_lists[j - n][0] != "G":
                        if rep_to_lists[j - n][0].startswith("S") == False and rep_to_lists[j - n][0].startswith(
                                "G") == False:
                            labelled_contigs[i].append("Undefined")
                            break
                        n += 1
                    if len(labelled_contigs[i]) < 5:
                        labelled_contigs[i].append(rep_to_lists[j - n][2])

                    n = 1
                    while rep_to_lists[j - n][0] != "F":
                        if rep_to_lists[j - n][0].startswith("S") == False and rep_to_lists[j - n][0].startswith(
                                "G") == False and rep_to_lists[j - n][0].startswith("F") == False:
                            labelled_contigs[i].append("Undefined")
                            break
                        n += 1
                    if len(labelled_contigs[i]) < 6:
                        labelled_contigs[i].append(rep_to_lists[j - n][2])
                ####################################################################################################

                elif actual_tax_level == "G":
                    labelled_contigs[i].append("---")
                    labelled_contigs[i].append(labelled_contigs[i][1])
                    n = 1
                    while rep_to_lists[j - n][0] != "F":
                        if rep_to_lists[j - n][0].startswith("S") == False and rep_to_lists[j - n][0].startswith(
                                "G") == False and rep_to_lists[j - n][0].startswith("F") == False:
                            labelled_contigs[i].append("Undefined")
                            break
                        n += 1
                    if len(labelled_contigs[i]) < 6:
                        labelled_contigs[i].append(rep_to_lists[j - n][2])
                ####################################################################################################
                elif actual_tax_level.startswith("F") and len(actual_tax_level) > 1:
                    # sono a livello sottospecie
                    labelled_contigs[i].append("---")
                    labelled_contigs[i].append("---")
                    n = 1
                    while rep_to_lists[j - n][0] != "F":
                        if rep_to_lists[j - n][0].startswith("S") == False and rep_to_lists[j - n][0].startswith(
                                "G") == False and rep_to_lists[j - n][0].startswith("F") == False:
                            labelled_contigs[i].append("Undefined")
                            break
                        n += 1
                    if len(labelled_contigs[i]) < 6:
                        labelled_contigs[i].append(rep_to_lists[j - n][2])
                ####################################################################################################
                elif actual_tax_level == "F":
                    labelled_contigs[i].append("---")
                    labelled_contigs[i].append("---")
                    labelled_contigs[i].append(labelled_contigs[i][1])
                ####################################################################################################
                else:
                    labelled_contigs[i].append("---")
                    labelled_contigs[i].append("---")
                    labelled_contigs[i].append("---")
