####################################################################################
# Script builds tibble of contrasts for the genera level and family level datasets #
####################################################################################

library(tidyverse)

## Begin get original pairs and PAML output for contrasts
# Not taking species richness and host variables b/c I recalculate them

genera_data <- read_csv("Raw_Sister_Tables/Papilionoidea_Genera_Contrasts.csv") %>%
    mutate(
        pair_id = as.character(Pair_no), dataset = "genera",
        taxa_left = Genus_first...2, taxa_right = Genus_last...3, # Note these columns were duplicated
        dN_left = dN_first, dN_right = dN_last,
        dS_left = dS_first, dS_right = dS_last,
        baseml_bl_left = blen_first, baseml_bl_right = blen_last
    ) %>%
    select(
        pair_id, dataset,
        taxa_left, taxa_right,
        dN_left, dN_right,
        dS_left, dS_right,
        baseml_bl_left, baseml_bl_right
    ) %>%
    pivot_longer(
        cols = ends_with("_left") | ends_with("_right"),
        names_to = c(".value", "side"),
        names_pattern = "(.*)_(left|right)"
    )

family_data <- read_csv("Raw_Sister_Tables/Lepidoptera_Families_Contrasts.csv") %>%
    mutate(
        pair_id = as.character(Pair_no), dataset = "family",
        taxa_left = Family_first, taxa_right = Family_last,
        dN_left = dN_first, dN_right = dN_last,
        dS_left = dS_first, dS_right = dS_last,
        baseml_bl_left = blen_first, baseml_bl_right = blen_last
    ) %>%
    select(
        pair_id, dataset,
        taxa_left, taxa_right,
        dN_left, dN_right,
        dS_left, dS_right,
        baseml_bl_left, baseml_bl_right
    ) %>%
    pivot_longer(
        cols = ends_with("_left") | ends_with("_right"),
        names_to = c(".value", "side"),
        names_pattern = "(.*)_(left|right)"
    )

major_lineage_data <- read_csv("Raw_Sister_Tables/Lepidoptera_MajorLineages_Contrasts.csv") %>%
    mutate(
        pair_id = as.character(Pair_no), dataset = "major-lineage",
        taxa_left = Lineage_first, taxa_right = Lineage_last
    ) %>%
    select(
        pair_id, dataset,
        taxa_left, taxa_right,
    ) %>%
    pivot_longer(
        cols = ends_with("_left") | ends_with("_right"),
        names_to = c(".value", "side"),
        names_pattern = "(.*)_(left|right)"
    )

pairs_long <- bind_rows(genera_data, family_data, major_lineage_data)
## End get original pairs and PAML output for contrasts


## Begin format  major lineage taxonomy
# To be joined into species count and host data below

# Convenience function to filter pseudo-duplicates which are an issue in the taxonomy
# Two rows with the same values except for NA in some places are pseudo-duplicates
filter_pseudo_duplicates <- function(df) {
    df <- distinct(df)
    if (dim(df)[1] <= 1) {
        return(df)
    }

    df_prime <- mutate(df, across(everything(), ~ replace_na(.x, "NA")))

    drop <- logical(nrow(df_prime))

    for (i in 1:(nrow(df_prime) - 1)) {
        if (drop[i]) next

        row <- unlist(df_prime[i, ])

        for (j in (i + 1):nrow(df_prime)) {
            if (drop[j]) next

            row_prime <- unlist(df_prime[j, ])

            if (all(row == row_prime)) {
                drop[j] <- TRUE
                next
            } else if (all(row != row_prime)) {
                next
            }

            unmatched <- which(!(row == row_prime))

            if (all(row[unmatched] == "NA")) {
                drop[i] <- TRUE
                break
            }
            if (all(row_prime[unmatched] == "NA")) drop[j] <- TRUE
        }
    }

    return(df[!drop, ])
}

taxonomy_major_lineage <- read_csv("Data/Lepi_MajorLineages/Lepi_MajorLineages_taxonomy.csv") %>% # Test
    filter(!if_all(everything(), ~ . == "")) %>%
    rename_with(tolower) %>%
    rowwise() %>%
    mutate(
        genus = str_to_title(str_trim(genus))
    ) %>%
    group_by(superfamily, family, genus) %>% # Filter duplicates within superfamily, family group. faster!
    group_modify(~ filter_pseudo_duplicates(.x)) %>%
    ungroup()

## End format major lineage taxonomy

## Begin wrangle host data
host_datasets <- c(
    "Data/Papi_Genera/Papi_Genera_HOSTS.csv",
    "Data/Lepi_MajorLineages/Lepi_MajorLineages_HOSTS.csv",
    "Data/Lepi_Families/Lepi_Families_HOSTS.csv"
)

# Note merging major-lineage taxonomy at end
host_data <- lapply(host_datasets, read_csv) %>%
    bind_rows() %>%
    filter(!if_all(everything(), ~ . == "")) %>%
    rename_with(tolower) %>%
    # Remove bad rows e.g. 'NA, Anarpia, CONNECTION_FAILED' (sic) in original
    filter(!(is.na(family) | is.na(genus) | is.na(species))) %>%
    mutate(
        host_name = str_to_title(str_trim(host_name)),
        host_family = str_to_title(str_trim(host_family))
    ) %>%
    mutate(
        genus = str_to_title(str_trim(genus)),
        family = str_to_title(str_trim(family)),
        host_genus = sub(" .*", "", host_name)
    ) %>%
    distinct() %>%
    group_by(species, host_genus) %>%
    filter(
        !(nchar(host_name) == nchar(host_genus) & n_distinct(host_name) > 1)
    ) %>%
    ungroup() %>%
    # Below should be many-to-one. But major-lineage taxonomy
    # has some species in two higher orders of the same level
    # I rely on joining with the chosen pairs to select the right one
    left_join(
        taxonomy_major_lineage,
        by = c("genus", "family"),
        relationship = "many-to-many" # Note below ^
    ) %>%
    # Here manually update a few problematic entries in Host data
    # TODO: Guard against same species in different higher orders?
    mutate(case_when(
        species == "Mnasilus allubitus" ~ "Papias allubitus",
        species == "Timoachis ruptifasciatus" ~ "Timochares ruptifasciata",
        species == "Cleosiris catamitus" ~ "Tetraganus catamitus",
        species == "Cleosiris lycaenoides" ~ "Tetraganus lycaenoides"
    ))


# Convenience function counts hosts at specified taxonomic level (group_var)
host_counts <- function(data, group_var) {
    data %>%
        group_by({{ group_var }}, species) %>%
        summarise(
            n_host_families = n_distinct(host_family),
            n_host_species = n_distinct(host_name)
        ) %>%
        mutate(is_generalist = ifelse(n_host_species > 1, 1, 0)) %>%
        ungroup() %>%
        group_by({{ group_var }}) %>%
        summarise(
            prop_generalist = sum(is_generalist) / n(),
            n_host_families = sum(n_host_families),
            n_host_species = sum(n_host_species)
        )
}

host_counts_genus <- host_counts(host_data, genus) %>%
    filter(genus != "Cosmopterigidae")

host_counts_family <- host_counts(host_data, family) %>% mutate(family = sub("^,", "", family))
host_counts_tribe <- host_counts(host_data, tribe)
host_counts_subfamily <- host_counts(host_data, subfamily)

combined_host_counts <- bind_rows(
    host_counts_genus, host_counts_family,
    host_counts_tribe, host_counts_subfamily
) %>%
    rowwise() %>%
    mutate(taxa = coalesce(genus, family, tribe, subfamily)) %>%
    select(taxa, prop_generalist, n_host_families, n_host_species)
## End wrangle host data

## Begin wrangle lepindex data
family_lepindex <- read_csv("Data/Papi_Genera/Lepi_Families_Lepindex_Raw.csv") %>%
    filter(!if_all(everything(), ~ . == "")) # Remove weird empty rows
papilionoidea_lepindex <- read_csv("Data/Papi_Genera/Papi_Genera_Lepindex_Raw.csv") %>%
    filter(!if_all(everything(), ~ . == ""))
major_lineage_lepindex <- read_csv("Data/Lepi_MajorLineages/Lepi_MajorLineages_Lepindex.csv") %>%
    select(Genus, Name_data) %>% # Avoid first number column
    filter(!if_all(everything(), ~ . == ""))

species_data <- bind_rows(
        family_lepindex, papilionoidea_lepindex, major_lineage_lepindex
    ) %>%
    rename_with(tolower) %>%
    rowwise() %>%
    filter(str_detect(name_data, "Valid Name")) %>%
    distinct(name_data, .keep_all = T) %>%
    mutate(species = word(name_data, 1)) %>%
    filter(genus != species) %>% # When genus matches species (avoid n() - 1 in original code)
    mutate(
        species = str_to_lower(species),
        genus = str_to_title(genus),
        family = str_to_title(family)
    ) %>%
    left_join(
        taxonomy_major_lineage,
        by = c("genus", "family"),
        relationship = "many-to-many" # See comment above on many-to-many
    ) %>%
    mutate(species = paste(genus, species)) %>%
    # Biblis genus case. All species should be hyperia
    mutate(species = ifelse(genus == "Biblis", "Biblis hyperia", species)) %>%
    # Drop rows with NA family and genus, species, name_data identical in another row
    group_by(genus, species, name_data) %>%
    filter(!(is.na(family) & n() > 1)) %>%
    ungroup() %>%
    # Record species with host in host_data
    left_join(
        host_data %>% select(species) %>% distinct() %>% mutate(host_record = 1),
        by = "species"
    ) %>%
    mutate(host_record = replace_na(host_record, 0))

species_counts <- function(data, group_var, level_name) {
    data %>%
        group_by({{ group_var }}) %>%
        summarise(
            n_species = n_distinct(species),
            n_host_record = sum(host_record)
        ) %>%
        mutate(level = level_name) %>%
        rename(taxa = {{ group_var }}) %>%
        drop_na()
}

species_counts_genus <- species_counts(species_data, genus, "genus")
species_counts_family <- species_counts(species_data, family, "family")
species_counts_tribe <- species_counts(species_data, tribe, "tribe")
species_counts_subfamily <- species_counts(species_data, subfamily, "subfamily")

combined_species_counts <- bind_rows(
        species_counts_genus, species_counts_family,
        species_counts_tribe, species_counts_subfamily
    ) %>%
    select(level, taxa, n_species, n_host_record)
## End wrangle lepindex data

## Begin associate recalculated host-related variables for contrasts, add weighting
pairs_long <- pairs_long %>%
    # TODO: Thorybes genus also in family host count?
        #   "Agaristinae" Should be taken care of by removing major-lineage duplicates
    left_join(combined_species_counts, by = "taxa") %>%
    left_join(combined_host_counts, by = "taxa") %>%
    mutate(
        weight_reciprocal = (n_species - n_host_record) / (n_host_record * (n_species - 1))
    )

# Sanity check
pairs_long %>% filter(n_host_record == 0 & n_host_species > 0)

# Calculate contrasts
contrasts <- pairs_long %>%
    group_by(dataset, pair_id) %>%
    mutate(taxa_size = n_species) %>%
    pivot_wider(
        names_from = side,
        values_from = -c(dataset, pair_id, side),
        names_glue = "{.value}_{side}",
        values_fill = NA
    ) %>%
    mutate(label = paste(taxa_left, "-", taxa_right)) %>%
    mutate(
        sign = ifelse(
            n_species_left - n_species_right != 0,
            sign(n_species_left - n_species_right),
            sample(c(1, -1), 1) # Randomise if tie (hence seed above)
        ),
        weight = 1 / (weight_reciprocal_left + weight_reciprocal_right)
    ) %>%
    select(-starts_with("weight_reciprocal")) %>%
    mutate(
        prop_generalist_left = asin(sqrt(prop_generalist_left)),
        prop_generalist_right = asin(sqrt(prop_generalist_right))
    ) %>%
    mutate(across(
        .cols = c(
            "dN_left", "dN_right", "dS_left", "dS_right",
            "n_host_species_left", "n_host_species_right",
            "n_species_left", "n_species_right",
            "n_host_families_left", "n_host_families_right",
            "baseml_bl_left", "baseml_bl_right"
        ),
        .fns = log
    )) %>%
    mutate(
        dN_contrast = sign * (dN_left - dN_right),
        dS_contrast = sign * (dS_left - dS_right),
        n_host_species_contrast = sign * (n_host_species_left - n_host_species_right),
        n_species_contrast = sign * (n_species_left - n_species_right),
        n_host_families_contrast = sign * (n_host_families_left - n_host_families_right),
        prop_generalist_contrast = sign * (prop_generalist_left - prop_generalist_right),
        baseml_bl_contrast = sign * (baseml_bl_left - baseml_bl_right)
    )

## End make contrast table

## Begin save contrasts and summary tables
write_csv(contrasts, file = "figures_and_output/contrasts.csv")
write_csv(combined_host_counts, file = "figures_and_output/combined_host_counts.csv")
write_csv(combined_species_counts, file = "figures_and_output/combined_species_counts.csv")

contrasts %>%
    group_by(dataset) %>%
    summarise(n = n()) %>%
    write_csv("figures_and_output/sample_sizes.csv")

## End
