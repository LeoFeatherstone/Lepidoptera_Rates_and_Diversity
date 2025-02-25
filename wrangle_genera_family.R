####################################################################################
# Script builds tibble of contrasts for the genera level and family level datasets #
####################################################################################

library(tidyverse)
library(phytools)
library(GGally)

# Convenience functions
source("R/baseml.extract.tree.R")

# Extracts left and right branch lengths from pair in PAML output
# Tree is a tree, measure refers to variable name. baseml_bl for baseml dN, dS for codeml
extract_paml_lengths <- function(tree, measure) {

  # Tree is a triplet, so sister pair always descends from the node with two descendants
  sister_parent_node <- which(node.depth(tree) == 2)
  pair <- tree$tip.label[getDescendants(tree, sister_parent_node)]
  left_taxa <- pair[1]
  right_taxa <- pair[2]

  terminal_edges <- tree$edge.length[tree$edge[, 2] <= Ntip(tree)]
  
  # Ensure we select from edges less than the number of tips
  left_blen <- terminal_edges[tree$tip.label == left_taxa]
  right_blen <- terminal_edges[tree$tip.label == right_taxa]
  
  result <- setNames(
    data.frame(
      left_taxa, left_blen,
      right_taxa, right_blen
    ),
    c("taxa_left", paste0(measure, "_left"), "taxa_right", paste0(measure, "_right"))
  )
  return(result)
}

# Extracts genus names from full taxa name
extract_genus <- function(name) {
    split <- unlist(str_split(name, "_"))
    split <- split[split != ""]
    taxa <- split[length(split) - 1]
    return(taxa)
}

# Extracts family names from full taxa name
extract_family <- function(name) {
    split <- unlist(str_split(name, "_"))
    split <- split[split != ""]
    taxa <- split[1]
    return(taxa)
}

# Extract branch length data from trees and bind to tibble
extract_and_bind <- function(trees, measure) {
    lengths <- lapply(trees, function(tree) extract_paml_lengths(tree, measure))
    df <- bind_rows(lengths, .id = "file") %>%
        mutate(
            dataset = ifelse(grepl("Papi_Genera", file), "genera", "family"),
            pair_id = str_extract(file, "\\d+")
        ) %>%
        select(-file)
    return(df)
}

# Prepare pecies diversity data. Handle family and genera species counts separately
family_lepindex <- read_csv("Data/Papi_Genera/Lepi_Families_Lepindex_Raw.csv") %>%
    filter(!if_all(everything(), ~ . == "")) # Remove weird empty rows
papilionoidea_lepindex <- read_csv("Data/Papi_Genera/Papi_Genera_Lepindex_Raw.csv") %>%
    filter(!if_all(everything(), ~ . == ""))

species_data <- bind_rows(family_lepindex, papilionoidea_lepindex) %>%
    rename_with(tolower) %>%
    filter(str_detect(name_data, "Valid Name")) %>%
    distinct(name_data, .keep_all = T)

species_counts_genus <- species_data %>%
    group_by(genus) %>% 
    summarise(n_species = n()) %>% # NB: distinct names filtered above
    mutate(genus = str_to_title(genus)) %>% # Lowercase to match pair_data$genus
    select(genus, n_species) %>%
    drop_na()

species_counts_family <- species_data %>%
    group_by(family) %>% 
    summarise(n_species = n()) %>% # NB: distinct names filtered above
    mutate(family = str_to_title(family)) %>% # Lowercase to match pair_data$family
    select(family, n_species) %>%
    drop_na()

combined_species_counts <- bind_rows(species_counts_genus, species_counts_family) %>%
    rowwise() %>%
    mutate(taxonomic_group = ifelse(is.na(genus), family, genus)) %>%
    select(taxonomic_group, n_species)

# Preparing host counts
host_datasets <- c(
        "Data/Papi_Genera/Papi_Genera_HOSTS.csv",
        "Data/Lepi_MajorLineages/Lepi_MajorLineages_HOSTS.csv",
        "Data/Lepi_Families/Lepi_Families_HOSTS.csv"
    )

# Key step below handles the case when genus-only samples are
# present alongside specific species. For example:
#    host_genus
#    Citrus
#    Citrus limonia
# Here, 'Citrus' counts as one host iff it is the only host of 
# a particular species from the genus citrus
host_data <- lapply(host_datasets, read_csv) %>%
    bind_rows() %>%
    filter(!if_all(everything(), ~ . == "")) %>% # Remove weird empty rows
    rename_with(tolower) %>%
    mutate(
        genus = str_to_title(str_trim(genus)),
        family = str_to_title(str_trim(family)),
        host_genus = sub(" .*", "", host_name)
    ) %>%
    distinct() %>%
    group_by(species, host_genus) %>%
    # Key step!
    filter(
        !(nchar(host_name) == nchar(host_genus) & n_distinct(host_name) > 1)
    ) %>%
    ungroup()

# Cout unique hosts species and families separately, then combine
host_counts_genus <- host_data %>%
    group_by(genus, species) %>%
    summarise(
        n_host_families = n_distinct(host_family),
        n_host_species = n_distinct(host_name)
    ) %>%
    mutate(is_generalist = ifelse(n_host_species > 1, 1, 0)) %>%
    ungroup() %>%
    group_by(genus) %>%
    summarise(
        prop_generalist = sum(is_generalist) / n(),
        n_host_families = sum(n_host_families),
        n_host_species = sum(n_host_species)
    ) %>%
    filter(genus != "Cosmopterigidae") # One off step. Possible error row in Papi hosts?

host_counts_family <- host_data %>%
    group_by(family, species) %>%
    summarise(
        n_host_families = n_distinct(host_family),
        n_host_species = n_distinct(host_name)
    ) %>%
    mutate(is_generalist = ifelse(n_host_species > 1, 1, 0)) %>%
    ungroup() %>%
    group_by(family) %>%
    summarise(
        prop_generalist = sum(is_generalist) / n(),
        n_host_families = sum(n_host_families),
        n_host_species = sum(n_host_species)
    )

combined_host_counts <- bind_rows(host_counts_genus, host_counts_family) %>%
    rowwise() %>%
    mutate(taxonomic_group = ifelse(is.na(genus), family, genus)) %>%
    select(taxonomic_group, prop_generalist, n_host_families, n_host_species)

# Parse PAML output  as lists of trees for each measure (baseml_bl, dN, dS)
# Handles genera and family level data simultaneously
baseml_trees_files <- list.files(
        path = c("PAML_outputs/Papi_Genera", "PAML_outputs/Lepi_Families"),
        pattern = "_baseml.txt", full.names = TRUE
    ) %>%
    str_sort(num = TRUE)

baseml_trees <- lapply(baseml_trees_files[-20], baseml.extract.tree) # NB. 20th family level codeml result is missing. TODO
names(baseml_trees) <- baseml_trees_files[-20]

codeml_trees_files <- list.files(
    path = c("PAML_outputs/Papi_Genera", "PAML_outputs/Lepi_Families"),
        pattern = "_codeml.txt", full.names = TRUE
    )  %>%
    str_sort(num = TRUE)

codeml_trees <- lapply(codeml_trees_files, codeml.extract.trees)
names(codeml_trees) <- codeml_trees_files
codeml_trees_dN <- lapply(codeml_trees, function(list) list$dN)
codeml_trees_dS <- lapply(codeml_trees, function(list) list$dS)

baseml_bl_df <- extract_and_bind(baseml_trees, "baseml_bl")
dN_df <- extract_and_bind(codeml_trees_dN, "dN")
dS_df <- extract_and_bind(codeml_trees_dS, "dS")

# Pair data ready to join species counts and host counts
pair_data_long <- baseml_bl_df %>%
    left_join(dS_df, by = c("pair_id", "dataset", "taxa_left", "taxa_right")) %>%
    left_join(dN_df, by = c("pair_id", "dataset", "taxa_left", "taxa_right")) %>%
    pivot_longer(
        cols = ends_with("_left") | ends_with("_right"),
        names_to = c(".value", "side"),
        names_pattern = "(.*)_(left|right)"
    ) %>% 
    rowwise() %>%
    mutate(
        taxonomic_group = ifelse(dataset == "family", extract_family(taxa), extract_genus(taxa))
    ) %>%
    mutate(taxonomic_group = ifelse(taxonomic_group == "NS0181", "Setabis", taxonomic_group)) # One-off

# Join species counts and host counts
#pair_data <- 

pair_data_wide <- pair_data_long %>%
    left_join(combined_species_counts, by = "taxonomic_group") %>%
    left_join(combined_host_counts, by = "taxonomic_group") %>%
    group_by(dataset, pair_id) %>%
    pivot_wider(
        names_from = side,
        values_from = -c(dataset, pair_id, side),
        names_glue = "{.value}_{side}",
        values_fill = NA
    )

# Get contrasts
contrasts <- pair_data_wide %>%
    group_by(dataset) %>%
    mutate(label = paste(taxonomic_group_left, "-", taxonomic_group_right)) %>% # TODO: For ggplotly later?
    mutate(across(
        .cols = -c(
            pair_id, taxa_left, taxa_right,
            taxonomic_group_left, taxonomic_group_right, label,
            prop_generalist_left, prop_generalist_right
        ),
        .fns = log
    )) %>%
    mutate(
        sign = sign(n_species_left - n_species_right)
    ) %>%
    mutate(
        prop_generalist_left = asin(sqrt(prop_generalist_left)),
        prop_generalist_right = asin(sqrt(prop_generalist_right))
    ) %>%
    mutate(
        dN_contrast = sign * (dN_left - dN_right),
        dS_contrast = sign * (dS_left - dS_right),
        n_host_species_contrast = sign * (n_host_species_left - n_host_species_right),
        n_species_contrast = sign * (n_species_left - n_species_right),
        n_host_families_contrast = sign * (n_host_families_left - n_host_families_right),
        baseml_bl_contrast = sign * (baseml_bl_left - baseml_bl_right),
        prop_generalist_contrast = sign * (prop_generalist_left - prop_generalist_right)
    ) %>%
    select(c(ends_with("contrast"), "label"))

# Save contrasts, host_data, species_data to .Rds file
saveRDS(
    list(
    contrasts = contrasts,
    host_data = host_data,
    species_data = species_data
    ),
    file = "contrasts_host_species_data.Rds"
)

## Plots (TODO: Must pull in major lineage data above for path analyses)
