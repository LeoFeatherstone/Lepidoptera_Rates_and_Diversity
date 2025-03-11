# Lepidopteran rates and diversity

This is a fork of [Lepidoptera_Rates_and_Diversity by amritchie](https://github.com/amritchie/Lepidoptera_Rates_and_Diversity). It takes the raw PAML output, recalculates host and species richness data, and performs analysis and plotting. I have included one new directory and two new scripts. All the original files and directories start with capitals, and my additions start with lowercase letters. 

## Reproducability
Users can reproduce the analysis by running each `_Run.R` file up to the points of writing the pair data (e.g., `Papilionoidea_Genera_Run.R`). Then run `wrangle_data.R`. From there, all plots and analyses can be reproduced in `plots_and_analysis.R`.

## Additions
- `/figures_and_output/`: Manuscript figures and tables for the manuscript.
- `wrangle_data.R`: Reads the original data for each pair in `Raw_Sister_Tables` and recalculates host and diversity variables.
- `plots_and_analysis.R`: Contains code to do path analyses and make plots.

## Original README

Below is description of the other directories and files from the original repository:

- `/R/`: Code to select sister pairs and produce alignments and contrasts for the three data sets used in this study. Each file can be sourced independently. If running from scratch, codeml and baseml results will need to be produced from alignments and pair trees and the file sourced again.

- `/Data/`: Initial sequence and phylogenetic data used for pair selection and production of pairwise alignments.

- `/Pair_alignments/`: Trimmed multiple sequence alignments for each sister pair plus outgroup. Can be input directly into PAML. There is no data for the MajorLineages data set because no substitution rates were calculated for this dataset.
- `/Pair_trees/`: Tree topologies for each pair, with all species names and sequence identifiers labelled on the tips. These are used as input topologies for PAML. There is no data for the MajorLineages dataset because PAML was not run and host data was taken from all species available for each tribe, subfamily, or family.

- `/PAML_outputs/`: Raw output files for baseml and codeml. Used to produce molecular branch length contrasts.

- `/R/`: Functions used in the main R code files.

- `/Raw_Sister_Tables/`: Edited CSV files with all output sister pair contrasts.

- `/Results/`: Final results as included in the manuscript.

