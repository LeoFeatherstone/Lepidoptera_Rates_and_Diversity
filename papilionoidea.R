# Papilionoidea data - wrangle, regress, plot

## TO DO:
# - codeml data - omega?
# - Count hosts per species
# - Possibly rerun PAML?

library(tidyverse)
library(phytools)
library(GGally)

# Convenience functions
source("R/baseml.extract.tree.R")
source("R/my.drop.tip.R")
source("R/ordered.ls.dir.R")
source("R/phylo.average.brlen.R")
source("R/welch.test.R")

get.cherries <- function(phy) {
    sapply(which(node.depth(phy) == 2), function(x) phy$tip.label[getDescendants(phy, x)], simplify = FALSE)
}

# Load sequence alignment and clean row names. TODO: Can we not make it slow by converting to matrix?
aln <- read.dna("Data/Papi_Genera/S_2c_Dataset_with_COI.phy", format = "sequential", as.matrix = TRUE)
rownames(aln) <- str_trim(rownames(aln))

# Load partition data
partition <- read.csv("Data/Papi_Genera/Papi_Genera_partitions.csv")

# Create a matrix indicating gaps in the alignment
gaps <- matrix(aln %in% c('a','c','g','t'), nrow = dim(aln)[1])
rownames(gaps) <- rownames(aln)

# Calculate coverage for each partition between first and last position
coverage <- mapply(
    function(first, last) gaps[, first:last] %*% rep(1, last - first + 1),
    partition$first, partition$last
)
colnames(coverage) <- partition$gid # gene ID
rownames(coverage) <- rownames(gaps)

# Select samples with complete coverage in genes of interest
genes <- colnames(coverage)[-which(colnames(coverage) %in% c("ArgKin", "COI-begin", "COI-end", "RpS2"))]
complete_taxa <- rowSums(coverage[, genes] > 0) == ncol(coverage[, genes])
coverage_complete <- coverage[complete_taxa, genes]
coverage_complete <- coverage_complete[rownames(coverage_complete) != "RWH_96_0878_Macrosoma_bahiata", ] # errant sample

# Load phylogenetic tree, filter to tips with complete coverage
# NB: Must get cherries, filter coverage, and THEN drop tips.
#     If you drop tips based on coverage first, you get artifical cherries!
tree <- read.nexus("Data/Papi_Genera/S_3b_core_analysis_median_ages.tre")
cherries <- get.cherries(tree)
pairs <- cherries[sapply(cherries, function(tips) all(tips %in% rownames(coverage_complete)))]

# Drop tip for outgroups later
# tree_complete_coverage <- drop.tip(
#     tree,
#     setdiff(tree$tip.label, rownames(coverage_complete))
# )


# Build data frame of paired taxa (pairs)
pairs <- lapply(pairs, function(pair) as.data.frame(t(pair)))
pairs <- bind_rows(pairs, .id = "id")
colnames(pairs) <- c("pair_id", "left_taxa", "right_taxa")

###### TODO, maybe build in PAML runs at this point ######

# Get omega, dN, dS, branch length from PAML output
baseml_trees_files <- str_sort(
    list.files("PAML_outputs/Papi_Genera", pattern = "_baseml.txt", full.names = TRUE),
    num = TRUE
)
baseml_trees <- lapply(baseml_trees_files, baseml.extract.tree)

codeml_trees_files <- str_sort(
    list.files("PAML_outputs/Papi_Genera", pattern = "_codeml.txt", full.names = TRUE),
    num = TRUE
)
codeml_trees <- lapply(codeml_trees_files, codeml.extract.trees)
codeml_trees_dN <- lapply(codeml_trees, function(list) list$dN)
codeml_trees_dS <- lapply(codeml_trees, function(list) list$dS)
codeml_trees_omega <- lapply(codeml_trees, function(list) list$w)

# Function to extract branch lengths for a given pair of taxa
# TODO: Add phylogenetic average case for family level case
extract_branch_lengths <- function(pair, tree, measure) {

  left_taxa <- pair$left_taxa
  right_taxa <- pair$right_taxa
  id <- pair$pair_id
  terminal_edges <- tree$edge.length[tree$edge[, 2] <= Ntip(tree)]
  
  # Ensure we select from edges less than the number of tips
  left_blen <- terminal_edges[tree$tip.label == left_taxa]
  right_blen <- terminal_edges[tree$tip.label == right_taxa]
  
  result <- setNames(
    data.frame(
      id,
      left_taxa, right_taxa,
      left_blen, right_blen
    ),
    c(
        "pair_id",
        "left_taxa", "right_taxa",
        paste0(measure, "_left"), paste0(measure, "_right")
    )
  )
}

# Get dN, dS, w, baseml data to vectors, then bind to df
# TODO: Omega has hashes in branch lengths??
baseml_bl <- list()
dN <- list()
dS <- list()
#omega <- list()
for (i in seq_along(pairs[, 1])) {
    baseml_bl[[i]] <- extract_branch_lengths(pairs[i, ], baseml_trees[[i]], "baseml_bl")
    dN[[i]] <- extract_branch_lengths(pairs[i, ], codeml_trees_dN[[i]], "dN")
    dS[[i]] <- extract_branch_lengths(pairs[i, ], codeml_trees_dS[[i]], "dS")
    # omega[[i]] <- extract_branch_lengths(pairs[i, ], codeml_trees_omega[[i]], "omega")
}
baseml_bl_df <- bind_rows(baseml_bl)
dN_df <- bind_rows(dN)
dS_df <- bind_rows(dS)

pair_data <- pairs %>%
    left_join(baseml_bl_df, by = c("pair_id", "left_taxa", "right_taxa")) %>%
    left_join(dN_df, by = c("pair_id", "left_taxa", "right_taxa")) %>%
    left_join(dS_df, by = c("pair_id", "left_taxa", "right_taxa"))

# Now add species counts to pair_data. Start with isolating genus

# Gets genus names from full taxa name
# TODO: Handle > tmp$name[100]
# [1] "NS0181_Setabis_"
extract_genus <- function(name) {
    split <- unlist(str_split(name, "_"))
    split <- split[split != ""]
    genus <- split[length(split) - 1]
    return(genus)
}

pair_data <- pair_data %>%
    rowwise() %>% # Ensure rowsie mutating. TODO: Why need this?
    mutate(
        genus_left = extract_genus(left_taxa),
        genus_right = extract_genus(right_taxa),
    ) %>%
    mutate(
        genus_left = ifelse(genus_left == "NS0181", "Setabis", genus_left),
        genus_right = ifelse(genus_right == "NS0181", "Setabis", genus_right)
    )

# Add species diversity data
family_lepindex <- read.csv("Data/Papi_Genera/Lepi_Families_Lepindex_Raw.csv", stringsAsFactors = F)
papilionoidea_lepindex <- read.csv("Data/Papi_Genera/Papi_Genera_Lepindex_Raw.csv", stringsAsFactors = F)

species_data <- as_tibble(bind_rows(family_lepindex, papilionoidea_lepindex)) %>%
    rename_with(tolower) %>%
    distinct(name_data, .keep_all = T) %>% 
    filter(str_detect(name_data, "Valid Name"))
    rename_with(tolower)

species_counts <- species_data %>%
    group_by(genus) %>% 
    summarise(n_species = n()) %>% # TODO: Check double counts like for hosts?
    mutate(genus = str_to_title(genus)) %>% # Lowercase to match pair_data$genus
    select(genus, n_species)

pair_data <- pair_data %>%
    pivot_longer(
        cols = starts_with("genus_"),
        values_to = "genus",
        names_to = "side"
    ) %>%
    left_join(species_counts, by = "genus") %>%
    mutate(side = gsub(side, pattern = "genus_", replacement = "")) %>%
    pivot_wider(
        names_from = side,
        values_from = c(n_species, genus),
        names_glue = "{.value}_{side}"
    )

# Now add host data. Wrangle to count host families and species
host_data <- read.csv("Data/Papi_Genera/Papi_Genera_HOSTS.csv")

# Key step below handles the case when genus-only samples are
# present alongside specific species. For example:
#    host_genus
#    Citrus
#    Citrus limonia
# Here, 'Citrus' counts as one host iff it is the only host of 
# a particular species from the genus citrus
host_data <- host_data %>%
    filter(!if_all(everything(), ~ . == "")) %>% # Remove weird empty rows
    rename_with(tolower) %>%
    mutate(
        genus = str_to_title(str_trim(genus)),
        host_genus = sub(" .*", "", host_name)
    ) %>%
    group_by(species, host_genus) %>%
    # Key step!
    filter(
        !(nchar(host_name) == nchar(host_genus) & n_distinct(host_name) > 1)
    ) %>%
    ungroup()

# Break point above touse `host_data` to look at num species with host data
host_counts <- host_data %>%
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
    )

pair_data <- pair_data %>%
    pivot_longer(cols = starts_with("genus_"), values_to = "genus", names_to = "side") %>%
    mutate(side = gsub(side, pattern = "genus_", replacement = "")) %>%
    left_join(host_counts, by = "genus") %>%
    pivot_wider(
        names_from = side,
        values_from = c(n_host_species, n_host_families, prop_generalist, genus),
        names_glue = "{.value}_{side}",
        values_fill = NA
    )

# Calculate contrasts and do pairs plot
contrasts <- pair_data %>%
    mutate(label = paste(genus_left, "-", genus_right)) %>%
    mutate(across(
        .cols = -c(
            pair_id, left_taxa, right_taxa,
            genus_left, genus_right, label,
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

# Pairs Plot
pairs_plot <- ggpairs(contrasts, columns = 1:7) +
    labs(
        title = "Log-transformed contrasts with pearson correlation",
        subtitle = "(prop_generalist is arcsin-sqrt transformed)"
    ) +
    theme(text = element_text(size = 11.5))

ggsave(plot = pairs_plot, filename = "genera_pairs.pdf", width = 12, height = 12, units = "in")

# n_species against dS
species_dS_plot <- contrasts %>%
    ggplot(aes(x = dS_contrast, y = n_species_contrast)) +
    geom_smooth(method= "lm", se = FALSE, col = alpha("black", 0.6)) +
    geom_point(fill = alpha("dodgerblue", 0.7), shape = 21) +
    ggpubr::stat_cor(method = "pearson", label.x = -0.4, label.y = 5) +
    labs(x = "dS contrast", y = "Number of species contrast", title = "Papilionoidea genera") +
    theme(
        text = element_text(size = 12)
    )
ggsave(
    plot = species_dS_plot, filename = "species_dS_genera.pdf",
    width = 4, height = 4, units = "in"
)

# Looking at recorded hosts for each species
matched_species_host <- species_data %>%
    mutate(
        species = word(name_data, 1),
        genus_species = paste(str_to_title(genus), species),
    ) %>%
    filter(!(genus == species)) %>% # Case when genus given as species at new genus
    select(c(family, genus, genus_species)) %>%
    rowwise() %>%
    left_join(
        host_data %>% group_by(species) %>% summarise(n_hosts = n_distinct(host_name)), 
        by = c("genus_species" = "species")
    ) %>%
    mutate(n_hosts = replace_na(n_hosts, 0))

# Shows that 99% of species don't have host data recorded! (98.9% technically...)
match_table <- matched_species_host %>%
    group_by(n_hosts) %>%
    summarise(
        num_species = n()
    ) %>%
    mutate(
        prop_species = num_species / sum(num_species)
    )

# Match table, finding most represented families among those with hose data
plot <- matched_species_host %>%
    group_by(family) %>%
    summarise(
        num_hosts_recorded = sum(n_hosts)
    ) %>%
    ungroup() %>%
    slice_max(num_hosts_recorded, n = 20) %>%
    drop_na() %>%
    arrange(desc(num_hosts_recorded)) %>%
    ggplot(aes(x = num_hosts_recorded, y = fct_reorder(family, num_hosts_recorded))) +
    geom_col(fill =alpha("dodgerblue", 0.75)) +
    labs(x = "Number of recorded hosts", y = "Family (top 20)") +
    theme_minimal() +
    theme(
        text = element_text(size = 14)
    )

ggsave(plot = plot, filename = "top_host_families.pdf", width = 6, height = 6, units = "in")

