##################################
# Plots and statistical analyses #
##################################

library(tidyverse)
library(GGally)
library(ggpubr)
library(lavaan)
library(igraph)

# Read the data
contrasts <- read_csv("contrasts.csv")
host_data <- read_csv("combined_host_counts.csv")
species_data <- read_csv("combined_species_counts.csv")

## Filter outliers
outlier_dN <- contrasts %>%
    filter(dataset == "genera") %>%
    slice_min(dN_contrast) %>%
    pull(label)

outlier_dS <- contrasts %>%
    filter(dataset == "family") %>%
    slice_max(dS_contrast) %>%
    pull(label)

## Remove outliers from contrasts
contrasts <- contrasts %>%
    filter(!(label %in% c(outlier_dN, outlier_dS)))

## Begin pairs plots
pairs_plot_pooled <- contrasts %>%
    rename_with(~ str_remove(., "_contrast")) %>%
    select(n_species, dN, dS, n_host_species, n_host_families, prop_generalist, dataset) %>%
    ggpairs(
        columns = 1:6,
        ggplot2::aes(colour = dataset),
        upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("smooth", alpha = 0.5, se = FALSE, lwd = 0.25)),
        diag = list(continuous = wrap("barDiag", alpha = 0.5, bins = 20))
    ) +
    labs(
        title = "Difference in log-transformed variables for all datasets. No outliers.",
        subtitle = "(prop_generalist is arcsin-sqrt transformed)"
    ) +
    scale_color_manual(values = c("red", "dodgerblue", "black")) +
    scale_fill_manual(values = c("red", "dodgerblue", "black")) +
    theme(text = element_text(size = 11.5))

ggsave(
    plot = pairs_plot_pooled,
    filename = "pairs_pooled_no_outliers.pdf",
    width = 12, height = 12, units = "in"
)
## End pairs plots

## Begin n_species against dS

# Fit linear regression models for each group with y-intercept fixed at 0
models <- contrasts %>%
    filter(dataset != "major-lineage") %>%
    group_by(dataset) %>%
    do(model = lm(n_species_contrast ~ dS_contrast, data = .))

# Extract coefficients, p-values, and correlations for each group
fit <- models %>%
    rowwise() %>%
    mutate(
        p_value = summary(model)$coefficients[1, 4],
        r_squared = summary(model)$r.squared
    ) %>%
    select(dataset, p_value, r_squared)

# Join model parameters into df
species_dS_plot <- contrasts %>%
    filter(dataset != "major-lineage") %>%
    left_join(fit, by = "dataset") %>%
    ggplot(aes(x = dS_contrast, y = n_species_contrast, fill = dataset, col = dataset)) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point(shape = 21, col = "black") +
    geom_label(
        aes(
            x = -0.1, y = 4.25,
            label = paste0("italic(p):", scales::scientific(p_value), "~italic(R^2):", scales::scientific(r_squared))
        ),
        col = "black",
        fill = "white",
        size = 3.5,
        parse = TRUE
    ) +
    scale_fill_manual(values = alpha(c("red", "dodgerblue"), 0.5)) +
    scale_color_manual(values = alpha(c("red", "dodgerblue"), 0.5)) +
    facet_wrap(
        ~dataset,
        scales = "free",
        labeller = as_labeller(c(family = "Lepidoptera families", genera = "Papilionoidea genera"))
    ) +
    labs(x = "Difference in log(dS)", y = "Difference in log(number of species)") +
    theme(
        legend.position = "none"
    )
ggsave(
    plot = species_dS_plot, filename = "species_vs_dS.pdf",
    width = 6, height = 3, units = "in"
)
## End n_species against dS


## Begin N species in taxa vs number of host records
host_vs_species_plot <- species_data %>%
    ggplot(aes(x = n_species, y = n_host_record)) +
    geom_point(shape = 21, fill = alpha("dodgerblue", 0.7)) +
    facet_wrap(
        ~level,
        scales = "free",
        labeller = as_labeller(c(
            family = "LepIndex families",
            genus = "LepIndex genera",
            subfamily = "LepIndex subfamilies",
            tribe = "LepIndex tribes"
        ))
    ) +
    labs(x = "Species in taxa", y = "Species with host record") +
    theme(
        text = element_text(size = 12)
    )

ggsave(
    plot = host_vs_species_plot, filename = "host_vs_species.pdf",
    width = 6, height = 3, units = "in"
)
## End N species in taxa vs number of host records

## Begin Welch and Waxman tests of ds
## Begin scatterplot for genera and family level data
welch_waxman <- contrasts %>%
    filter(dataset %in% c("genera", "family")) %>%
    mutate(
        abs_dS_scaled = abs(dS_left - dS_right) / sqrt(abs(baseml_bl_left - baseml_bl_right)),
        sqrt_abs_bl = sqrt(abs(baseml_bl_left - baseml_bl_right))
    ) %>%
    ggplot(aes(x = sqrt_abs_bl, y = abs_dS_scaled, fill = dataset, col = dataset)) +
    geom_point(shape = 21, col = "black") +
    scale_fill_manual(values = alpha(c("red", "dodgerblue"), 0.5)) +
    scale_color_manual(values = alpha(c("red", "dodgerblue"), 0.5)) +
    labs(x = expression(sqrt("|baseml contrast|")), y = expression(frac("|dS contrast|", sqrt("|baseml contrast|")))) +
    theme(
        text = element_text(size = 12),
        legend.position = "none"
    ) +
    facet_wrap(~dataset, scales = "free", labeller = as_labeller(c(family = "Lepidoptera families", genera = "Papilionoidea genera")))
ggsave(
    plot = welch_waxman, filename = "welch-waxman.pdf",
    width = 6, height = 3, units = "in"
)
## End scatterplot for genera and family level data
## End Welch and Waxman

# Define the path models
# NB If convergence issue, remove dN as a direct effect on N species
path_model_genera_family <- "
  # Direct effects
  # Note dN as an effect below. Removed due to convergence issues
  n_species ~ n_host_species + n_host_families + prop_generalist + dS

  # Indirect effect
  dN ~ n_host_species + n_host_families + prop_generalist + b * dS

  # Covariances
  n_host_families ~~ c_host_fam_species * n_host_species
  prop_generalist ~~ c_gen_host_species * n_host_species

  # Constraints
  c_host_fam_species > 0
  c_gen_host_species > 0
  b > 0

"

path_model_major_lineage <- "
  # Direct effects
  n_species ~ n_host_species + n_host_families + prop_generalist

  # Covariances
  n_host_families ~~ c_host_fam_species * n_host_species
  prop_generalist ~~ c_gen_host_species * n_host_species

  # Constraints
  c_host_fam_species > 0
  c_gen_host_species > 0

"

# Fit the model for each group
genera_model <- contrasts %>%
    rename_with(~ str_remove(., "_contrast")) %>%
    filter(dataset == "genera") %>%
    mutate(weight = ifelse(is.na(weight), 0, weight)) %>%
    as.data.frame() %>%
    sem(data = ., model = path_model_genera_family, sampling.weights = "weight")

# TODO: Weight for major lineage data
major_lineage_model <- contrasts %>%
    rename_with(~ str_remove(., "_contrast")) %>%
    filter(dataset == "major-lineage") %>%
    mutate(weight = ifelse(is.na(weight), 0, weight)) %>%
    as.data.frame() %>%
    sem(data = ., model = path_model_major_lineage, sampling.weights = "weight")

family_model <- contrasts %>%
    rename_with(~ str_remove(., "_contrast")) %>%
    filter(dataset == "family") %>%
    mutate(weight = ifelse(is.na(weight), 0, weight)) %>%
    as.data.frame() %>%
    sem(data = ., path_model_genera_family, sampling.weights = "weight")

fit_models <- list(genera_model, major_lineage_model, family_model)

# Convenience function to make network plots for each model with igraph
plot_model <- function(model, title) {
    set.seed(1234) # To ensure same layout. 1000 looks good

    # Extract parameter estimates
    params <- parameterestimates(model) %>%
        filter(op == "~" | op == "~~") %>%
        mutate(significant = ifelse(pvalue < 0.05, 1, 0)) %>%
        mutate(
            lhs = recode(lhs, 
                n_species = "Species", 
                n_host_species = "Host\nspecies", 
                n_host_families = "Host\nfamilies", 
                prop_generalist = "Proportion\ngeneralist", 
                dS = "dS", 
                dN = "dN"
            ),
            rhs = recode(rhs, 
                n_species = "Species", 
                n_host_species = "Host\nspecies", 
                n_host_families = "Host\nfamilies", 
                prop_generalist = "Proportion\ngeneralist", 
                dS = "dS", 
                dN = "dN"
            )
        )

    # Create edge list and filter self-loops
    edges <- params %>%
        select(lhs, rhs, est, significant, op) %>%
        rename(from = rhs, to = lhs, weight = est) %>%
        filter(from != to)

    # Create graph object
    g <- graph_from_data_frame(edges, directed = TRUE)
    layout = layout_in_circle(g)

    # Define edge colors based on the effect
    edge_colors <- case_when(
        edges$significant == 1 & edges$weight < 0 ~ "dodgerblue",
        edges$significant == 1 & edges$weight > 0 ~ "orange",
        TRUE ~ alpha("black", 0.4)
    )

    # Define edge shapes based on the operation
    edge_shapes <- ifelse(edges$op == "~~", 3, 1)

    # Plot graph with fixed layout
    plot(
        g,
        layout = layout,
        edge.label = round(E(g)$weight, 2),
        edge.color = edge_colors,
        edge.label.color = "black",
        edge.arrow.mode = edge_shapes,

        vertex.label.color = "black",
        vertex.label.font = 2,
        vertex.shape = "circle",
        vertex.size = 45,
        vertex.color = "white",
        vertex.label.family = "Helvetica",
        main = title
    )
}

pdf("path-analyses.pdf", width = 12, height = 12)

par(mfrow = c(2, 2))
plot_model(genera_model, "Papilonoidea genera")
plot_model(major_lineage_model, "Lepidoptera major lineages")
plot_model(family_model, "Lepidoptera families")

# Add legend for weight color
plot.new()
legend(
    "center",
    legend = c("Positive significant", "Negative significant", "Non-significant"),
    col = c("orange", "dodgerblue", alpha("black", 0.4)), lty = 1, cex = 1.2
)

dev.off()

# Fit the model for each group without weights
genera_model_no_weights <- contrasts %>%
    rename_with(~ str_remove(., "_contrast")) %>%
    filter(dataset == "genera") %>%
    as.data.frame() %>%
    sem(data = ., model = path_model_genera_family)

major_lineage_model_no_weights <- contrasts %>%
    rename_with(~ str_remove(., "_contrast")) %>%
    filter(dataset == "major-lineage") %>%
    as.data.frame() %>%
    sem(data = ., model = path_model_major_lineage)

family_model_no_weights <- contrasts %>%
    rename_with(~ str_remove(., "_contrast")) %>%
    filter(dataset == "family") %>%
    as.data.frame() %>%
    sem(data = ., model = path_model_genera_family)

fit_models_no_weights <- list(genera_model_no_weights, major_lineage_model_no_weights, family_model_no_weights)

# Plot models without weights
pdf("path-analyses-no-weights.pdf", width = 12, height = 12)

par(mfrow = c(2, 2))
plot_model(genera_model_no_weights, "Papilonoidea genera (no weights)")
plot_model(major_lineage_model_no_weights, "Lepidoptera major lineages (no weights)")
plot_model(family_model_no_weights, "Lepidoptera families (no weights)")

# Add legend for weight color
plot.new()
legend(
    "center",
    legend = c("Positive significant", "Negative significant", "Non-significant"),
    col = c("orange", "dodgerblue", alpha("black", 0.4)), lty = 1, cex = 1.2
)

dev.off()

# Fit the model for each group with baseml_bl instead of dN and dS
path_model_genera_family_bl <- "
    # Direct effects
    n_species ~ n_host_species + n_host_families + prop_generalist + baseml_bl

    # Indirect effect
    baseml_bl ~ n_host_species + n_host_families + prop_generalist

    # Covariances
    n_host_families ~~ c_host_fam_species * n_host_species
    prop_generalist ~~ c_gen_host_species * n_host_species

    # Constraints
    c_host_fam_species > 0
    c_gen_host_species > 0
"

path_model_major_lineage_bl <- "
    # Direct effects
    n_species ~ n_host_species + n_host_families + prop_generalist

    # Covariances
    n_host_families ~~ c_host_fam_species * n_host_species
    prop_generalist ~~ c_gen_host_species * n_host_species

    # Constraints
    c_host_fam_species > 0
    c_gen_host_species > 0
"

# Fit the model for each group with weights
genera_model_bl <- contrasts %>%
        rename_with(~ str_remove(., "_contrast")) %>%
        filter(dataset == "genera") %>%
        mutate(weight = ifelse(is.na(weight), 0, weight)) %>%
        as.data.frame() %>%
        sem(data = ., model = path_model_genera_family_bl, sampling.weights = "weight")

major_lineage_model_bl <- contrasts %>%
        rename_with(~ str_remove(., "_contrast")) %>%
        filter(dataset == "major-lineage") %>%
        mutate(weight = ifelse(is.na(weight), 0, weight)) %>%
        as.data.frame() %>%
        sem(data = ., model = path_model_major_lineage_bl, sampling.weights = "weight")

family_model_bl <- contrasts %>%
        rename_with(~ str_remove(., "_contrast")) %>%
        filter(dataset == "family") %>%
        mutate(weight = ifelse(is.na(weight), 0, weight)) %>%
        as.data.frame() %>%
        sem(data = ., model = path_model_genera_family_bl, sampling.weights = "weight")

fit_models_bl <- list(genera_model_bl, major_lineage_model_bl, family_model_bl)

# Plot models with weights
pdf("path-analyses-bl.pdf", width = 12, height = 12)

par(mfrow = c(2, 2))
plot_model(genera_model_bl, "Papilonoidea genera (baseml_bl)")
plot_model(major_lineage_model_bl, "Lepidoptera major lineages (baseml_bl)")
plot_model(family_model_bl, "Lepidoptera families (baseml_bl)")

# Add legend for weight color
plot.new()
legend(
        "center",
        legend = c("Positive significant", "Negative significant", "Non-significant"),
        col = c("orange", "dodgerblue", alpha("black", 0.4)), lty = 1, cex = 1.2
)

dev.off()

# Fit the model for each group without weights
genera_model_no_weights_bl <- contrasts %>%
        rename_with(~ str_remove(., "_contrast")) %>%
        filter(dataset == "genera") %>%
        as.data.frame() %>%
        sem(data = ., model = path_model_genera_family_bl)

major_lineage_model_no_weights_bl <- contrasts %>%
        rename_with(~ str_remove(., "_contrast")) %>%
        filter(dataset == "major-lineage") %>%
        as.data.frame() %>%
        sem(data = ., model = path_model_major_lineage_bl)

family_model_no_weights_bl <- contrasts %>%
        rename_with(~ str_remove(., "_contrast")) %>%
        filter(dataset == "family") %>%
        as.data.frame() %>%
        sem(data = ., model = path_model_genera_family_bl)

fit_models_no_weights_bl <- list(genera_model_no_weights_bl, major_lineage_model_no_weights_bl, family_model_no_weights_bl)

# Plot models without weights
pdf("path-analyses-no-weights-bl.pdf", width = 12, height = 12)

par(mfrow = c(2, 2))
plot_model(genera_model_no_weights_bl, "Papilonoidea genera (no weights, baseml_bl)")
plot_model(major_lineage_model_no_weights_bl, "Lepidoptera major lineages (no weights, baseml_bl)")
plot_model(family_model_no_weights_bl, "Lepidoptera families (no weights, baseml_bl)")

# Add legend for weight color
plot.new()
legend(
        "center",
        legend = c("Positive significant", "Negative significant", "Non-significant"),
        col = c("orange", "dodgerblue", alpha("black", 0.4)), lty = 1, cex = 1.2
)

dev.off()
