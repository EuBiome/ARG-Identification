import argparse

ap = argparse.ArgumentParser()

ap.add_argument("--famiglia_path", required=True, help="Insert the path of the list of pathogens family")
ap.add_argument("--genere_path", required=True, help="Insert the path of the list of pathogens genre")
ap.add_argument("--specie_path", required=True, help="Insert the path of the list of pathogens species")
ap.add_argument("--microbial_report_path", required=True, help="Insert the path of kraken2 report")
ap.add_argument("--taxa_patogeni_path", required=True, help="Insert the path where you want to store the whole list of pathogens taxa")
args = vars(ap.parse_args())

famiglia_path = args["famiglia_path"]
genere_path = args["genere_path"]
specie_path = args["specie_path"]
microbial_report_path = args["microbial_report_path"]
taxa_patogeni_path = args["taxa_patogeni_path"]

with open(famiglia_path, mode="r") as famiglia_file, open(genere_path, mode="r") as genere_file,\
        open(specie_path, mode="r") as specie_file, open(microbial_report_path, mode="r") as microbial_report_file,\
        open(taxa_patogeni_path, mode="w") as taxa_patogeni_file:

    # Trasformo il file del report in una lista di liste:
    # Ogni sottolista è così strutturata: ex. ['G', '10494', 'Lymphocystivirus']
    # Perciò: (livello tassonomico, id_tassonomico, nomenclatura_tassonomica)
    report_matrix = []
    for line in microbial_report_file:
        report_to_list = [line.split()[3], line.split()[4], line.split("\t")[5].lstrip()[:-1]]
        report_matrix.append(report_to_list)

    famiglia_taxa_patogeni = []
    genere_taxa_patogeni = []
    specie_taxa_patogeni = []
    tutti_taxa_patogeni = []

    for famiglia in famiglia_file:
        famiglia_taxa_patogeni.append(famiglia.split("\t")[2][:-1])

    for genere in genere_file:
        genere_taxa_patogeni.append(genere.split("\t")[2][:-1])

    for specie in specie_file:
        specie_taxa_patogeni.append(specie.split("\t")[2][:-1])

    dentro_la_famiglia = False

    for taxid in famiglia_taxa_patogeni:
        for i in range(len(report_matrix)):
            if taxid == report_matrix[i][1]:
                tutti_taxa_patogeni.append(taxid)
                dentro_la_famiglia = True
            elif dentro_la_famiglia == True:
                if report_matrix[i][0].startswith("R") or report_matrix[i][0].startswith("D") or report_matrix[i][0].startswith("P") or report_matrix[i][0].startswith("O") or report_matrix[i][0].startswith("C") or report_matrix[i][0] == "F":
                    dentro_la_famiglia = False
                else:
                    tutti_taxa_patogeni.append(report_matrix[i][1])

    print(len(tutti_taxa_patogeni))

    dentro_il_genere = False

    for taxid in genere_taxa_patogeni:
        for i in range(len(report_matrix)):
            if taxid == report_matrix[i][1]:
                tutti_taxa_patogeni.append(taxid)
                dentro_il_genere = True
            elif dentro_il_genere == True:
                if report_matrix[i][0].startswith("R") or report_matrix[i][0].startswith("D") or report_matrix[i][
                    0].startswith("P") or report_matrix[i][0].startswith("O") or report_matrix[i][0].startswith("C") or \
                        report_matrix[i][0].startswith("F") or report_matrix[i][0] == "G":
                    dentro_il_genere = False
                else:
                    tutti_taxa_patogeni.append(report_matrix[i][1])

    print(len(tutti_taxa_patogeni))
    dentro_la_specie = False

    for taxid in specie_taxa_patogeni:
        for i in range(len(report_matrix)):
            if taxid == report_matrix[i][1]:
                tutti_taxa_patogeni.append(taxid)
                dentro_la_specie = True
            elif dentro_la_specie == True:
                if report_matrix[i][0].startswith("R") or report_matrix[i][0].startswith("D") or report_matrix[i][
                    0].startswith("P") or report_matrix[i][0].startswith("O") or report_matrix[i][0].startswith("C") or \
                        report_matrix[i][0].startswith("F") or report_matrix[i][0].startswith("G") or report_matrix[i][0] == "S":
                    dentro_la_specie = False
                else:
                    tutti_taxa_patogeni.append(report_matrix[i][1])

    tutti_taxa_patogeni = list(set(tutti_taxa_patogeni))
    
    for taxa in tutti_taxa_patogeni:
        taxa_patogeni_file.write(taxa + "\n")
