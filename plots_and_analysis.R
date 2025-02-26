##################################
# Plots and statistical analyses #
##################################

library(tidyverse)
library(GGally)
library(ggpubr)
# Read the data
contrasts_host_species_data <- readRDS("contrasts_host_species_data.Rds")

contrasts <- contrasts_host_species_data$contrasts
host_data <- contrasts_host_species_data$host_data
species_data <- contrasts_host_species_data$species_data

# Pairs plots
pairs_plot_genus <- contrasts %>%
    filter(dataset == "genera") %>%
    ungroup() %>%
    select(-c("dataset", "label", "weight", "n_species_left", "n_species_right", starts_with("n_host_record"))) %>%
    ggpairs() +
    labs(
        title = "Papilonoidea genera: log-transformed contrasts with pearson correlation",
        subtitle = "(prop_generalist is arcsin-sqrt transformed)"
    ) +
    theme(text = element_text(size = 11.5))

pairs_plot_major_lineage <- contrasts %>%
    filter(dataset == "major-lineage") %>%
    ungroup() %>%
    select(-c("dataset", "label", "dN_contrast", "dS_contrast", "baseml_bl_contrast", "weight", "n_species_left", "n_species_right", starts_with("n_host_record"))) %>%
    ggpairs() +
    labs(
        title = "Lepidoptera major-lineage: log-transformed contrasts with pearson correlation",
        subtitle = "(prop_generalist is arcsin-sqrt transformed)"
    ) +
    theme(text = element_text(size = 11.5))

pairs_plot_family <- contrasts %>%
    filter(dataset == "family") %>%
    ungroup() %>%
    select(-c("dataset", "label", "weight", "n_species_left", "n_species_right", starts_with("n_host_record"))) %>%
    ggpairs() +
    labs(
        title = "Lepidoptera family: log-transformed contrasts with pearson correlation",
        subtitle = "(prop_generalist is arcsin-sqrt transformed)"
    ) +
    theme(text = element_text(size = 11.5))

ggsave(plot = pairs_plot_genus, filename = "pairs_genera.pdf", width = 12, height = 12, units = "in")
ggsave(plot = pairs_plot_major_lineage, filename = "pairs_major-lineage.pdf", width = 12, height = 12, units = "in")
ggsave(plot = pairs_plot_family, filename = "pairs_family.pdf", width = 12, height = 12, units = "in")

# n_species against dS

# Fit linear regression models for each group with y-intercept fixed at 0
models <- contrasts %>%
  filter(dataset != "major-lineage") %>%
  group_by(dataset) %>%
  do(model = lm(n_species_contrast ~ 0 + dS_contrast, data = .))

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
    ggplot(aes(x = dS_contrast, y = n_species_contrast)) +
    geom_smooth(formula = y ~ x + 0, method= "lm", se = FALSE, col = alpha("black", 0.6)) +
    geom_point(fill = alpha("dodgerblue", 0.7), shape = 21) +
    geom_label(
        aes(
            x = -0.4, y = 4.25, 
            label = paste0("italic(p):", round(p_value, 2), "~italic(R^2):", round(r_squared, 2))
        ), 
        size = 3.5,
        parse = TRUE
    ) +
    facet_wrap(
        ~dataset, scales = "free",
        labeller = as_labeller(c(family = "Lepidoptera families", genera = "Papilionoidea genera"))
    ) +
    labs(x = "dS contrast", y = "Number of species contrast") +
    theme(
        text = element_text(size = 12)
    )
ggsave(
    plot = species_dS_plot, filename = "species_vs_dS.pdf",
    width = 6, height = 3, units = "in"
)

# N species in taxa vs number of host records
host_vs_species_plot <- species_data %>%
    ggplot(aes(x = n_species, y = n_host_record)) +
    geom_point(shape = 21, fill = alpha("dodgerblue", 0.7)) +
    facet_wrap(
        ~level, scales = "free",
        labeller = as_labeller(c(family = "LepIndex families", genus = "LepIndex genera"))
    ) +
    labs(x = "Number of species in family/genera", y = "Number with host record") +
    theme(
        text = element_text(size = 12)
    )

ggsave(
    plot = host_vs_species_plot, filename = "host_vs_species.pdf",
    width = 6, height = 3, units = "in"
)

# Prepare lavaan plot code
library(lavaan)
library(semPlot)
library(gridExtra)

# Define the path model
path_model <- '
    n_species_contrast ~ dS_contrast + dN_contrast + baseml_bl_contrast
    dS_contrast ~ dN_contrast + baseml_bl_contrast
    dN_contrast ~ baseml_bl_contrast
'

# Fit the model for each group
fit_models <- contrasts %>%
    group_by(dataset) %>%
    do(fit = sem(path_model, data = .))

# Plot the path models
plot_list <- fit_models %>%
    rowwise() %>%
    mutate(plot = list(semPaths(fit, "std", layout = "circle", title = dataset))) %>%
    pull(plot)

# Arrange plots side by side
grid.arrange(grobs = plot_list, ncol = 2)

# Save the combined plot
ggsave(
    filename = "path_models_combined.pdf",
    plot = grid.arrange(grobs = plot_list, ncol = 2),
    width = 16, height = 8, units = "in"
)