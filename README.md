# Lepidopteran rates and diversity

This is a fork of [Lepidoptera_Rates_and_Diversity by amritchie](https://github.com/amritchie/Lepidoptera_Rates_and_Diversity). It takes the raw PAML output, recalculates host and species richness data, and performs analysis and plotting. I have included one new directory and two new scripts. All the original files and directories start with capitals, and my additions start with lowercase letters. 

## Reproducibility
Users can reproduce the analysis by running each `_Run.R` file up to the points of writing the pair data (e.g., `Papilionoidea_Genera_Run.R`). Then run `wrangle_data.R`. From there, all plots and analyses can be reproduced in `plots_and_analysis.R`.

## Additions
- `/figures_and_output/`: Manuscript figures and tables for the manuscript.
- `wrangle_data.R`: Reads the original data for each pair in `Raw_Sister_Tables` and recalculates host and diversity variables.
- `plots_and_analysis.R`: Contains code to do path analyses and make plots.

## How does `wrangle_data.R` work?
The key part of this addition to the original work is `wrangle_data.R`, which reads the original data for each pair in `Raw_Sister_Tables` and recalculates host and diversity variables. Here is a step-by-step breakdown of how it works:

### 1 Parse original contrasts 
- Parse original contrast data for genera, family, and major lineage (L10-50)
     - These are `Raw_Sister_Tables/Papilionoidea_Genera_Contrasts.csv`, `Raw_Sister_Tables/Papilionoidea_Family_Contrasts.csv`, and `Raw_Sister_Tables/MajorLineages_Contrasts.csv`.
     - Only the PAML outputs are retained (dN, dS, and branch lengths)
- Parse major lineage taxonomy file (`Raw_Sister_Tables/Lepidoptera_MajorLineages_Contrasts.csv`) which is used to find the taxonomic level for pairs in the major lineages dataset (Family, Subfamily, or Tribe, L52-66).
    - Pseudo-duplicate rows are removed (L74-123). For example:

    | Family      | Tribe      | Superfamily |
    |-------------|------------|-------------|
    | Nymphalidae | Satyrini   | Papilionoidea|
    | Nymphalidae | NA         | Papilionoidea|
    
    Would be reduced to:

    | Family      | Tribe      | Superfamily |
    |-------------|------------|-------------|
    | Nymphalidae | Satyrini   | Papilionoidea|

### 2 Parse host data
- Parse and merge host datasets:
    - `Raw_Sister_Tables/Papilionoidea_Genera_Hosts.csv`
    - `Raw_Sister_Tables/Papilionoidea_Family_Hosts.csv`
    - `Raw_Sister_Tables/MajorLineages_Hosts.csv`
- Filter out erroneous rows. Here's an example of an erroneous row from the original:
    'NA, Anarpia, CONNECTION_FAILED' (*sic*)
- Within each combination of Lepi species and host species, filter duplicates where one row gives only the genus of the host and the other gives the full species. This is not done if there is only genus information for a particular host. Here is an example:

| Lepidoptera species | Host species |
|---------------------|--------------|
| Papilio nobilis | Warburgia |
|Papilio nobilis | Warburgia salutaris |

is condensed to:

| Lepidoptera species | Host species |
|---------------------|--------------|
|Papilio nobilis | Warburgia salutaris |

- Species and genera are renamed to match Lepindex taxonomy using the following rules:
    - `"Mnasilus allubitus"` → `"Papias allubitus"`
    - `"Timochares ruptifasciatus"` → `"Timochares ruptifasciata"`
    - `"Cleosiris catamitus"` → `"Tetragonus catamitus"`
    - `"Cleosiris lycaenoides"` → `"Tetragonus lycaenoides"`
    - Genus `"Mnasilus"` → `"Papias"`
    - Genus `"Cleosiris"` → `"Tetragonus"`

- Host diversity variables are now calculated via the following procedure. This is done by the `host_counts()` function in lines 183-218:
    1. If counting at the genera level, group the dataset by individual species and count the number of unique host genera and host families
    2. Mark which species are generalists (i.e., have more than one unique host genus or family)
    3. Then group at the level of interest (genera, family, or major lineage) and sum the counts of unique host genera and families. Also calculate the proportion of generalists at this level

### 3 Parse species richness data
- Parse and merge individual species:
    - `Raw_Sister_Tables/Papilionoidea_Genera_Species_Richness.csv`
    - `Raw_Sister_Tables/Papilionoidea_Family_Species_Richness.csv`
    - `Raw_Sister_Tables/MajorLineages_Species_Richness.csv`

- Filter down to only rows with "Valid name" status
- Manually add Biblis hyperia based on manual check, similar to above
- Filter out rows with repeated and incomplete informatio below family level (L252-253). For example rows such as the following would be removed:
    | Family      | Genus         |
    |-------------|---------------|
    | NA          | BUCCULATRIX   |
    | NA          | BUCCULATRIX   |
    | NA          | BUCCULATRIX   |

- Now merge host data into this dataset matching on lepidopteran species name. In doing so, add a variable called `host_record` (L257) that is a binary indicator of whether a given species has a host recorded in the host dataset. This is later used for weighting in path analysis

- Count the species richness in each genera, family, or major lineage by taking the sum of the number of distinct species at each level. This is done by the `species_counts()` function (L262-284).

### 4 Merge all data and make contrast table (L320-434)
- Merge the host and species richness data with the orignal contrasts for PAML output
- Once again manually adjust host and generalist measurements for a small subset of data that are not properly matched between the Lepindex and HOSTS datasets.

In a small number of cases, host and generalist variables are manually set for taxa that do not match cleanly between datasets. This is done using conditional logic in the code, for example:

```r
prop_generalist = case_when(
    taxa == "Doidae" ~ 1.0,
    taxa == "Sematuridae" ~ 2/3,
    TRUE ~ prop_generalist
),
n_host_families = case_when(
    taxa == "Doidae" ~ 1,
    taxa == "Epicopeiidae" ~ 7,
    taxa == "Sematuridae" ~ 9,
    TRUE ~ n_host_families
),
n_host_species = case_when(
    taxa == "Doidae" ~ 1,
    taxa == "Sematuridae" ~ 11,
    TRUE ~ n_host_species
),
n_host_record = case_when(
    taxa == "Doidae" ~ 1,
    taxa == "Sematuridae" ~ 3,
    TRUE ~ n_host_record
)
```

- Finally, log transform all variables of interest and calculate contrasts. If there is a tie between values for any variables **before** log transformation, then the sign of the contrast is randomised. Hence, rerunning this script repeatedly and with different seeds may results in a slightly different contrast table.


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

