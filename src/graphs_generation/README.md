By running the main script is possible to generate the following type of graphs starting from the kraken2 output and report and from the output of rgi. The graph relates the taxids of the bacteria with the antibiotics they are resistant to. It's possible to choose the graph at species or at genre level.
<p align="center">
  <img src="/imgs/Analysis_figure_g_g.png" width="800" title="Analysis_figure" alt="Analysis_figure">
</p>

## Example Usage

```
python3 /path/to/main.py --level_of_interest s --rgi_output /path/to/rgi.txt --kraken2_output /path/to/k2_results --rgi_analysis=/path/where/to/store/the/output --kraken2_report /path/to/k2_report
```
