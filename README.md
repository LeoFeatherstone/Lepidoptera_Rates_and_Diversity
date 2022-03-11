# Lepidoptera_Rates_and_Diversity
Code, data and results for relating Lepidopteran diversity, molecular rates and host plant variables.

This repository contains:

Top level - R code to select sister pairs and produce alignments and contrasts for the three data sets use din this study. Each file can be sourced independently. If running from scratch, codeml and baseml results will need to be produced from alignments and pair trees and the file sourced again.

Data - Initial sequence and phylogenetic data used for pair selection and production of pairwise alignments.

Pair_alignments - Trimmed multiple sequence alignments for each sister pair plus outgroup. Can be input directly into PAML. There is no data for the MajorLineages data set because no substitution rates were calculated for this dataset.

Pair_trees - Tree topologies for each pair, with all species names and sequence identifiers labelled on the tips. These are used as input toplogies for PAML. There is no data for the MajorLineages dataset because PAML was not run and host data was taken from all species available for each tribe, subfamily or family.

PAML_outputs - raw output files for baseml and codeml. Used to produce molecular branch length contrasts.

R - functions used in the main R code files.

Raw_Sister_Tables - edited CSV files with all output sister pair contrasts.

Results - Final results as included in the manuscript.
