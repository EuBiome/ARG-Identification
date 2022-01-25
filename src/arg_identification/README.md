By running the main inside this folder we execute the workflow shown in figure. The result is the production of the 4 lists, each one representing a different risk level of the detected arg genes.

<p align="center">
  <img src="/imgs/args_workflow.png" width="400" title="Args workflow" alt="Args workflow">
</p>

The build_pathogens_list.py script is usefull to produce the list of all pathogens taxid starting from 3 lists of pathogens at family, genre and species level and from the report_file of kraken2. If you already have this list you can directly skip to the main script.
