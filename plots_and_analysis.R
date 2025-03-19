##################################
# Plots and statistical analyses #
##################################

library(tidyverse)
library(GGally)
library(ggpubr)
library(lavaan)
library(igraph)

# Read the data
contrasts <- read_csv("figures_and_output/contrasts.csv")
host_data <- read_csv("figures_and_output/combined_host_counts.csv")
species_data <- read_csv("figures_and_output/combined_species_counts.csv")

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
    filename = "figures_and_output/pairs_pooled_no_outliers.pdf",
    width = 12, height = 12, units = "in", dpi = 300
)
ggsave(
    plot = pairs_plot_pooled,
    filename = "figures_and_output/pairs_pooled_no_outliers.png",
    width = 12, height = 12, units = "in", dpi = 600
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
    plot = species_dS_plot, filename = "figures_and_output/species_vs_dS.pdf",
    width = 6, height = 3, units = "in", dpi = 300
)
ggsave(
    plot = species_dS_plot, filename = "figures_and_output/species_vs_dS.png",
    width = 6, height = 3, units = "in", dpi = 300
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
    plot = host_vs_species_plot, filename = "figures_and_output/host_vs_species.pdf",
    width = 6, height = 3, units = "in", dpi = 300
)
ggsave(
    plot = host_vs_species_plot, filename = "figures_and_output/host_vs_species.png",
    width = 6, height = 3, units = "in", dpi = 300
)
## End N species in taxa vs number of host records

## Begin Welch and Waxman tests of ds
## Begin scatterplot for genera and family level data
## /sqrt(sum of branch length/2)
welch_waxman <- contrasts %>%
    filter(dataset %in% c("genera", "family")) %>%
    mutate(
        dS_left = exp(dS_left),
        dS_right = exp(dS_right),
        baseml_bl_left = exp(baseml_bl_left),
        baseml_bl_right = exp(baseml_bl_right)
    ) %>%
    mutate(
        abs_dS_scaled = abs(dS_left - dS_right) / sqrt((baseml_bl_left + baseml_bl_right) / 2),
        sqrt_sum_bl = sqrt((baseml_bl_left + baseml_bl_right) / 2)
    ) %>%
    ggplot(aes(x = sqrt_sum_bl, y = abs_dS_scaled, fill = dataset, col = dataset)) +
    geom_point(shape = 21, col = "black") +
    scale_fill_manual(values = alpha(c("red", "dodgerblue"), 0.5)) +
    scale_color_manual(values = alpha(c("red", "dodgerblue"), 0.5)) +
    labs(x = expression(sqrt("Sum of branch lengths")), y = expression(frac("|dS contrast|", sqrt("Sum of branch lengths")))) +
    theme(
        text = element_text(size = 12),
        legend.position = "none"
    ) +
    facet_wrap(~dataset, scales = "free", labeller = as_labeller(c(family = "Lepidoptera families", genera = "Papilionoidea genera")))
ggsave(
    plot = welch_waxman, filename = "figures_and_output/welch-waxman.pdf",
    width = 6, height = 3, units = "in", dpi = 300
)
ggsave(
    plot = welch_waxman, filename = "figures_and_output/welch-waxman.png",
    width = 6, height = 3, units = "in", dpi = 300
)
## End scatterplot for genera and family level data
## End Welch and Waxman

# Define the path models
# NB If convergence issue, remove dN as a direct effect on N species

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
                dN = "dN",
                baseml_bl = "Branch\nlength"
            ),
            rhs = recode(rhs,
                n_species = "Species",
                n_host_species = "Host\nspecies",
                n_host_families = "Host\nfamilies",
                prop_generalist = "Proportion\ngeneralist",
                dS = "dS",
                dN = "dN",
                baseml_bl = "Branch\nlength"
            )
        )

    # Create edge list and filter self-loops
    edges <- params %>%
        select(lhs, rhs, est, significant, op) %>%
        rename(from = rhs, to = lhs, weight = est) %>%
        filter(from != to)

    # Create graph object
    g <- graph_from_data_frame(edges, directed = TRUE)
    layout <- layout_in_circle(g)

    # Define edge colors based on the effect
    edge_colors <- case_when(
        edges$significant == 1 & edges$weight < 0 ~ "dodgerblue",
        edges$significant == 1 & edges$weight > 0 ~ "red",
        .default = "grey"
    )

    # Define edge shapes based on the operation
    edge_shapes <- ifelse(edges$op == "~~", 3, 1)

    # Plot graph with fixed layout
    plot(
        g,
        layout = layout,
        margin = c(0.1, 0.1, 0.1, 0.1),
        edge.label = round(E(g)$weight, 2),
        edge.arrow.mode = edge_shapes,
        edge.arrow.size = 0.5,
        edge.color = edge_colors,
        edge.label.color = "black",
        edge.label.cex 	= .75,
        edge.curved = -0.2,

        vertex.label.color = "black",
        vertex.shape = "circle",
        vertex.size = 60,
        vertex.label.cex = 0.75,
        vertex.label.font = 1,
        vertex.color = "white",
        vertex.label.family = "Helvetica",
        main = title
    )
}

path_model_genera_family <- "
  # Direct effects
  # Note dN as an effect below. Removed due to convergence issues
  n_species ~ n_host_species + n_host_families + prop_generalist

  # Covariances
  n_host_families ~~ n_host_species
  prop_generalist ~~ n_host_species
  n_host_families ~~ baseml_bl
  prop_generalist ~~ baseml_bl
  n_host_species ~~ baseml_bl
  n_species ~~ baseml_bl
"

path_model_major_lineage <- "
  # Direct effects
  # Note dN as an effect below. Removed due to convergence issues
  n_species ~ n_host_species + n_host_families + prop_generalist

  # Covariances
  n_host_families ~~ n_host_species
  prop_generalist ~~ n_host_species
"

# Fit the model for each group
genera_model <- contrasts %>%
    rename_with(~ str_remove(., "_contrast")) %>%
    filter(dataset == "genera") %>%
    mutate(weight = ifelse(is.na(weight), 0, weight)) %>%
    as.data.frame() %>%
    sem(data = ., model = path_model_genera_family, sampling.weights = "weight")

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

pdf("figures_and_output/path-analyses-final.pdf", width = 6, height = 6, useDingbats = FALSE)

    par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
    plot_model(genera_model, "Papilonoidea genera")
    plot_model(major_lineage_model, "Lepidoptera major lineages")
    plot_model(family_model, "Lepidoptera families")

    # Add legend for weight color
    plot.new()
    legend(
        "center",
        legend = c("Positive significant", "Negative significant", "Non-significant"),
        col = c("red", "dodgerblue", "grey"), lty = 1, cex = 1.2
    )

dev.off()

png("figures_and_output/path-analyses-final.png", width = 6, height = 6, units = "in", res = 600)

    par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
    plot_model(genera_model, "Papilonoidea genera")
    plot_model(major_lineage_model, "Lepidoptera major lineages")
    plot_model(family_model, "Lepidoptera families")

    # Add legend for weight color
    plot.new()
    legend(
        "center",
        legend = c("Positive significant", "Negative significant", "Non-significant"),
        col = c("red", "dodgerblue", "grey"), lty = 1, cex = 1.2
    )

dev.off()

# Fit family with no weights for comparison
family_model_no_weights <- contrasts %>%
    rename_with(~ str_remove(., "_contrast")) %>%
    filter(dataset == "family") %>%
    as.data.frame() %>%
    sem(data = ., model = path_model_genera_family)

pdf("figures_and_output/family-path-unweighted-final.pdf", width = 7, height = 4)
    par(mfrow = c(1, 2), mar = c(1, 1, 2.5, 1))
    plot_model(family_model_no_weights, "Lepidoptera families\n(unweighted)")
    plot.new()
    legend(
        "center",
        legend = c("Positive significant", "Negative significant", "Non-significant"),
        col = c("red", "dodgerblue", "grey"), lty = 1, cex = 0.8
    )
dev.off()

png("figures_and_output/family-path-unweighted-final.png",  width = 7, height = 4, units = "in", res = 600)
    par(mfrow = c(1, 2), mar = c(1, 1, 2.5, 1))
    plot_model(family_model_no_weights, "Lepidoptera families\n(unweighted)")
    plot.new()
    legend(
        "center",
        legend = c("Positive significant", "Negative significant", "Non-significant"),
        col = c("red", "dodgerblue", "grey"), lty = 1, cex = 0.8
    )
dev.off()


# Save parameter estimates table
# Extract parameter estimates from each model and write to a table
parameter_estimates <- lapply(
    c(fit_models, family_model_no_weights),
    parameterestimates
)
names(parameter_estimates) <- c("genera", "major_lineage", "family", "family-unweighted")

# Combine all parameter estimates into a single data frame and rename 'op' to 'direction'
parameter_estimates_df <- bind_rows(parameter_estimates, .id = "model") %>%
    rename(direction = op) %>%
    mutate(direction = recode(direction, `~` = "direct effect", `~~` = "covariance"))

# Write the parameter estimates to a CSV file
write_csv(parameter_estimates_df, "figures_and_output/path_estimates.csv")
