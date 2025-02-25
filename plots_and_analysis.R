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
    select(-c("dataset", "label")) %>%
    ggpairs() +
    labs(
        title = "Log-transformed contrasts with pearson correlation",
        subtitle = "(prop_generalist is arcsin-sqrt transformed)"
    ) +
    theme(text = element_text(size = 11.5))

pairs_plot_family <- contrasts %>%
    filter(dataset == "family") %>%
    ungroup() %>%
    select(-c("dataset", "label")) %>%
    ggpairs() +
    labs(
        title = "Log-transformed contrasts with pearson correlation",
        subtitle = "(prop_generalist is arcsin-sqrt transformed)"
    ) +
    theme(text = element_text(size = 11.5))

ggsave(plot = pairs_plot_genus, filename = "pairs_genera.pdf", width = 12, height = 12, units = "in")
ggsave(plot = pairs_plot_family, filename = "pairs_family.pdf", width = 12, height = 12, units = "in")

# n_species against dS

# Fit linear regression models for each group with y-intercept fixed at 0
models <- contrasts %>%
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
