By running the main inside this folder we execute the workflow shown in figure. The result is the production of the 4 lists, each one representing a different risk level of the detected arg genes.

<p align="center">
  <img src="/imgs/args_workflow.png" width="400" title="Args workflow" alt="Args workflow">
</p>

The build_pathogens_list.py script is usefull to produce the list of all pathogens taxid starting from 3 lists of pathogens at family, genre and species level and from the report_file of kraken2. If you already have this list you can directly skip to the main script.

## Example Usage

```
python3 /path/to/main.py --plasmid_classification_path /path/to/plasmid_classification_file --microbial_classification_path /path/to/microbial_classification_file --microbial_report_path /path/to/microbial_report_file --rgi_path /path/to/rgi_output_file --genes_list /path/to/AR_genes_list --pathogens_list /path/to/pathogens_list
```
```
python3 /path/to/build_pathogens_list.py --family_path /path/to/pathogens_families --genre_path /path/to/pathogens_genres --species_path /path/to/pathogens_species --microbial_classification_path /path/to/microbial_classification_file --microbial_report_path /path/to/microbial_report_file --pathogens_list /path/where/to/store/the/output
```
