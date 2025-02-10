library(tidyverse)
library(ggpubr)
library(cowplot)
library(lavaan)
library(lavaanPlot)
library(DiagrammeRsvg)
library(rsvg)
library(GGally)
library(igraph)
library(ggraph)

sister_table_files <- list.files(
    path = "Raw_Sister_Tables", pattern = "*.csv",
    full.names = TRUE
)
sister_table_list <- lapply(
    sister_table_files, function(x) read.csv(x, row.names = NULL)
)

names(sister_table_list) <- c(
    "Family", "Major lineage", "Genera"
)
sister_table_list %>%
    bind_rows(.id = "id") %>%
    colnames()

set.seed(1234)
data1 <- sister_table_list[[1]] %>%
    mutate_if(is.numeric, log) %>%
    mutate(sign = ifelse(sign(N_spp_first - N_spp_last) != 0, sign(N_spp_first - N_spp_last), sample(c(1, -1), 1, 0.5))) %>%
    mutate(
        dN = sign * (dN_first - dN_last),
        dS = sign * (dS_first - dS_last),
        n_spp = sign * (N_spp_first - N_spp_last),
        host_pgen = sign * (asin(sqrt(exp(Hosts_pgen_first))) - asin(sqrt(exp(Hosts_pgen_last)))),
        host_mean = sign * (Hosts_mean_first - Hosts_mean_last),
        host_families = sign * (Host_families_first - Host_families_last),
        host_species = sign * (Host_species_first - Host_species_last),
        host_pd = sign * (Hosts_pd_first - Hosts_pd_last)
    ) %>%
    select(
        host_pd, dN, dS, host_pgen, host_mean, host_families, host_species, n_spp
    )

data2 <- sister_table_list[[2]] %>%
    mutate_if(is.numeric, log) %>%
    mutate(sign = ifelse(sign(N_spp_first - N_spp_last) != 0, sign(N_spp_first - N_spp_last), sample(c(1, -1), 1, 0.5))) %>%
    mutate(
        #dN = sign * (dN_first - dN_last),
        #dS = sign * (dS_first - dS_last),
        n_spp = sign * (N_spp_first - N_spp_last),
        host_pgen = sign * (asin(sqrt(exp(Hosts_pgen_first))) - asin(sqrt(exp(Hosts_pgen_last)))),
        host_mean = sign * (Hosts_mean_first - Hosts_mean_last),
        host_families = sign * (Host_families_first - Host_families_last),
        host_species = sign * (Host_species_first - Host_species_last),
        host_pd = sign * (Hosts_pd_first - Hosts_pd_last)
    ) %>%
    # select(
    #     host_pd, dN, dS, host_pgen, host_mean, host_families, host_species, n_spp
    # )
    select(
        host_pd, host_pgen, host_mean, host_families, host_species, n_spp
    )

data3 <- sister_table_list[[3]] %>%
    mutate_if(is.numeric, log) %>%
    mutate(sign = ifelse(sign(N_spp_first - N_spp_last) != 0, sign(N_spp_first - N_spp_last), sample(c(1, -1), 1, 0.5))) %>%
    mutate(
        dN = sign * (dN_first - dN_last),
        dS = sign * (dS_first - dS_last),
        n_spp = sign * (N_spp_first - N_spp_last),
        host_pgen = sign * (asin(sqrt(exp(Hosts_pgen_first))) - asin(sqrt(exp(Hosts_pgen_last)))),
        host_mean = sign * (Hosts_mean_first - Hosts_mean_last),
        host_families = sign * (Host_families_first - Host_families_last),
        host_species = sign * (Host_species_first - Host_species_last),
        host_pd = sign * (Hosts_pd_first - Hosts_pd_last)
    ) %>%
    select(
        host_pd, dN, dS, host_pgen, host_mean, host_families, host_species, n_spp
    )

p1 <- ggpairs(data1, columns = 1:8)
p2 <- ggpairs(data2, columns = 1:6)
p3 <- ggpairs(data3, columns = 1:8)

p3

cowplot::plot_grid(p1, p2, p3, nrow = 1)

sister_table_list[[1]] %>%
    select(starts_with("Hosts_pd"))




## Example plots for rationale document

# Mockup random data
set.seed(123)
random_data <- data.frame(
    clade_size = rnorm(300, mean = 0, sd = 10),
    branch_length = rnorm(300, mean = 0, sd = 10),
    dN = rnorm(300, mean = 0, sd = 10),
    dS = rnorm(300, mean = 0, sd = 10),
    host_families = rnorm(300, mean = 0, sd = 10),
    host_species = rnorm(300, mean = 0, sd = 10),
    prop_generalist = runif(300, min = 0, max = 1),
    mean_hosts = rnorm(300, mean = 0, sd = 10),
    host_pd = rnorm(300, mean = 0, sd = 10),
    grouping = rep(c("Papilionoidea genera", "Lepidoptera major lineages", "Lepidoptera families"), each = 100)
)

p1 <- random_data %>%
    ggplot(aes(y = clade_size, x = branch_length)) +
    geom_point(shape = 21, fill = alpha("dodgerblue", 0.7)) +
    stat_cor(method = "kendall", cor.coef.name = "tau") +
    facet_wrap(~grouping)

p2 <- random_data %>%
    ggplot(aes(y = clade_size, x = dN)) +
    geom_point(shape = 21, fill = alpha("dodgerblue", 0.7)) +
    stat_cor(method = "kendall", cor.coef.name = "tau") +
    facet_wrap(~grouping)

p3 <- random_data %>%
    ggplot(aes(y = clade_size, x = dS)) +
    geom_point(shape = 21, fill = alpha("dodgerblue", 0.7)) +
    stat_cor(method = "kendall", cor.coef.name = "tau") +
    facet_wrap(~grouping)

p <- plot_grid(p1, p2, p3, align = "v", ncol = 1)

ggsave(p, file = "example_fig1.pdf", dpi = 600)

p4 <- random_data %>%
    ggplot(aes(y = clade_size, x = host_families)) +
    geom_point(shape = 21, fill = alpha("dodgerblue", 0.7)) +
    stat_cor(method = "kendall", cor.coef.name = "tau") +
    facet_wrap(~grouping)

p5 <- random_data %>%
    ggplot(aes(y = clade_size, x = host_species)) +
    geom_point(shape = 21, fill = alpha("dodgerblue", 0.7)) +
    stat_cor(method = "kendall", cor.coef.name = "tau") +
    facet_wrap(~grouping)

p6 <- random_data %>%
    ggplot(aes(y = clade_size, x = prop_generalist)) +
    geom_point(shape = 21, fill = alpha("dodgerblue", 0.7)) +
    stat_cor(method = "kendall", cor.coef.name = "tau") +
    facet_wrap(~grouping)

p7 <- random_data %>%
    ggplot(aes(y = clade_size, x = mean_hosts)) +
    geom_point(shape = 21, fill = alpha("dodgerblue", 0.7)) +
    stat_cor(method = "kendall", cor.coef.name = "tau") +
    facet_wrap(~grouping)

p8 <- random_data %>%
    ggplot(aes(y = clade_size, x = host_pd)) +
    geom_point(shape = 21, fill = alpha("dodgerblue", 0.7)) +
    stat_cor(method = "kendall", cor.coef.name = "tau") +
    facet_wrap(~grouping)

fig2 <- plot_grid(p4, p5, p6, p7, p8, align = "v", ncol = 1)

ggsave(fig2, file = "example_fig2.pdf", dpi = 600)

# Mock path analysis

# Define the model
model <- "
    # direct effects
    clade_size ~  dN + host_species + mean_hosts

    # indirect effects
    dN ~ dS
    host_species ~~ mean_hosts
"
fit <- sem(model, data = random_data)

# Plot the path diagram
path_plot <- lavaanPlot(
    model = fit,
    coefs = FALSE
)

# Save the combined plot
save_png(path_plot, "path_analysis_plot.png")


# Graph examples

# Define the lavaan model
model <- "
    # direct effects
    clade_size ~  dN + host_species + mean_hosts

    # indirect effects
    dN ~ dS
    host_species ~~ mean_hosts
"

# Manually parse the model to extract relationships
edges <- data.frame(
  from = c("dN", "host_species", "mean_hosts", "dS", "host_species"),
  to = c("clade_size", "clade_size", "clade_size", "dN", "mean_hosts"),
  relationship = c("direct", "direct", "direct", "direct", "correlation")
)

# Create a graph object
graph <- graph_from_data_frame(edges, directed = TRUE)

# Customize edge types for correlation
E(graph)$type <- edges$relationship

# Plot the graph using ggraph
ggraph(graph, layout = "stress") +
  geom_edge_link(aes(linetype = type), edge_width = 1, edge_color = "gray") +
  geom_node_point(size = 5, color = "steelblue") +
  geom_node_text(aes(label = name), vjust = 1.5, size = 5) +
  theme_void()

# Load required libraries
library(dplyr)

# Define all nodes in the graph
all_nodes <- c("dN", "host_species", "mean_hosts", "clade_size")

# Define the scenarios with their specific edges
scenarios <- list(
  "Incompatibility due to isolation" = data.frame(
    from = c("dN"),
    to = c("clade_size"),
    relationship = c("direct")
  ),
  "Escape and radiate" = data.frame(
    from = c("host_species"),
    to = c("clade_size"),
    relationship = c("direct")
  ),
  "Oscillate" = data.frame(
    from = c("mean_hosts", "host_species"),
    to = c("clade_size", "clade_size"),
    relationship = c("positive", "positive")
  ),
  "Musical chairs" = data.frame(
    from = c("mean_hosts", "host_species"),
    to = c("clade_size", "clade_size"),
    relationship = c("negative", "positive")
  )
)

# Create a complete set of possible connections
all_connections <- expand.grid(
  from = all_nodes,
  to = all_nodes,
  stringsAsFactors = FALSE
) %>%
  filter(from != to)  # Remove self-loops

# Add each scenario to the connections data
edges <- bind_rows(lapply(names(scenarios), function(scenario) {
  scenario_edges <- scenarios[[scenario]]
  
  # Join with all_connections to ensure every possible edge is included
  all_connections %>%
    left_join(scenario_edges, by = c("from", "to")) %>%
    mutate(scenario = scenario)  # Add scenario column
}))

edges %>%
  group_by(scenario) %>%
  do({
    graph <- graph_from_data_frame(., vertices = data.frame(name = all_nodes), directed = TRUE)
    ggraph(graph, layout = "stress") +
      geom_edge_link(aes(linetype = relationship), edge_width = 1, edge_color = "gray") +
      geom_node_point(size = 5, color = "steelblue") +
      geom_node_text(aes(label = name), vjust = 1.5, size = 5) +
      theme_void() +
      ggtitle(unique(.$scenario))
  }) %>%
  plot_grid(plotlist = ., labels = "AUTO", ncol = 2)



pivot_longer()