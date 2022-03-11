library(ape)
library(phytools)
library(tidyverse)
library(fs)
library(phylotate)
library(reshape2)

source("R\\baseml.extract.tree.R")
source("R\\my.drop.tip.R")
source("R\\ordered.ls.dir.R")
source("R\\phylo.average.brlen.R")
source("R\\welch.test.R")

lepi.annot.trees.rot <- sapply(dir_ls("Data\\Lepi_MajorLineages\\trees", regex=".*tre$"), read_annotated, simplify = F)

#lepi.annot.trees.rot <- sapply(lepi.annot.trees, function (x) if (any(str_detect(x$node.comment[!is.na(x$node.comment)], "rotate=true"))) {multi2di(rotateNodes(x, which(str_detect(x$node.comment, "rotate=true"))))} else {multi2di(x)}, simplify = F)

lepi_trees_names <- sapply(dir_ls("Data\\Lepi_MajorLineages\\trees", regex=".*txt$"), function (x) read.csv(x, header=F, stringsAsFactors=F), simplify = F)

lepi_taxonomy <- read.csv("Data\\Lepi_MajorLineages\\Lepi_MajorLineages_taxonomy.csv", stringsAsFactors = F)

lepi_tree_lineage_ts <- list()
for (i in seq_along(lepi.annot.trees.rot)) 
  {
  lepi.annot.trees.rot[[i]]$tip.label <- lepi_trees_names[[i]][,2][match(lepi.annot.trees.rot[[i]]$tip.label, lepi_trees_names[[i]][,1])]
  write.tree(lepi.annot.trees.rot[[i]], paste0("Pair_trees\\Lepi_MajorLineages\\Full_named_trees\\lepi_majorlineages_named_", i, ".tre"))

  tmp.gen <- toupper(str_split(lepi.annot.trees.rot[[i]]$tip.label, "_", simplify=T)[,1])
  tmp.lin <- sapply(tmp.gen, function (y) if (y %in% lepi_taxonomy$Genus) (lepi_taxonomy %>% filter(Genus == y))[1,] else tibble(Superfamily="UNKNOWN", Family="UNKNOWN", Subfamily="UNKNOWN", Tribe="UNKNOWN", Genus=y), simplify = F)

  lepi_tree_lineage_ts[[i]] <- bind_rows(tmp.lin)
  }

lepi_tree_tribes_tokeep <- list()
lepi_tree_subfam_tokeep <- list()
lepi_tree_family_tokeep <- list()
for (i in seq_along(lepi_tree_lineage_ts)) 
  {
    lepi_tree_lineage_ts[[i]]$Lineage  <- apply(lepi_tree_lineage_ts[[i]], 1, function (x) paste(x[!is.na(x)][1:length(x[!is.na(x)])-1], collapse="_"))
    lepi_tree_lineage_ts[[i]]$Tip <- lepi.annot.trees.rot[[i]]$tip.label
    lepi_tree_lineage_ts[[i]]$TribeIsMono <- sapply(lepi_tree_lineage_ts[[i]]$Tribe, function (x) if (is.na(x)) NA else is.monophyletic(lepi.annot.trees.rot[[i]], lepi_tree_lineage_ts[[i]] %>% filter(Tribe == x) %>% pull(Tip)))
    lepi_tree_lineage_ts[[i]]$SubfamIsMono <- sapply(lepi_tree_lineage_ts[[i]]$Subfamily, function (x) if (is.na(x)) NA else is.monophyletic(lepi.annot.trees.rot[[i]], lepi_tree_lineage_ts[[i]] %>% filter(Subfamily == x) %>% pull(Tip)))
    lepi_tree_lineage_ts[[i]]$FamIsMono <- sapply(lepi_tree_lineage_ts[[i]]$Family, function (x) if (is.na(x)) NA else is.monophyletic(lepi.annot.trees.rot[[i]], lepi_tree_lineage_ts[[i]] %>% filter(Family == x) %>% pull(Tip)))
    
    lepi_tree_tribes_tokeep[[i]] <- c(lepi_tree_lineage_ts[[i]] %>% filter(!is.na(Tribe) & TribeIsMono) %>% group_by(Tribe) %>% summarise(Tip = first(Tip)) %>% pull(Tip),lepi_tree_lineage_ts[[i]] %>% filter(!is.na(Tribe) & !TribeIsMono) %>% pull(Tip))
    lepi_tree_subfam_tokeep[[i]] <- c(lepi_tree_lineage_ts[[i]] %>% filter(!is.na(Subfamily) & SubfamIsMono) %>% group_by(Subfamily) %>% summarise(Tip = first(Tip)) %>% pull(Tip),lepi_tree_lineage_ts[[i]] %>% filter(!is.na(Subfamily) & !SubfamIsMono) %>% pull(Tip))
    lepi_tree_family_tokeep[[i]] <- c(lepi_tree_lineage_ts[[i]] %>% filter(!is.na(Family) & FamIsMono) %>% group_by(Family) %>% summarise(Tip = first(Tip)) %>% pull(Tip),lepi_tree_lineage_ts[[i]] %>% filter(!is.na(Family) & !FamIsMono) %>% pull(Tip))

  }

lepi.annot.trees.tribes <- sapply(seq_along(lepi_tree_lineage_ts), function (x) if (length(lepi_tree_tribes_tokeep[[x]])==0) NA else keep.tip(lepi.annot.trees.rot[[x]], lepi_tree_tribes_tokeep[[x]]), simplify = F)
lepi.annot.trees.subfam <- sapply(seq_along(lepi_tree_lineage_ts), function (x) keep.tip(lepi.annot.trees.rot[[x]], lepi_tree_subfam_tokeep[[x]]), simplify = F)
lepi.annot.trees.family <- sapply(seq_along(lepi_tree_lineage_ts), function (x) keep.tip(lepi.annot.trees.rot[[x]], lepi_tree_family_tokeep[[x]]), simplify = F)

#for (i in seq_along(lepi.annot.trees.tribes)) 
 # {
  # lepi.annot.trees.subfam[[i]]$tip.label <- as.character(lepi.annot.trees.subfam[[i]]$tip.label)
    #lepi.annot.trees.family[[i]]$tip.label <- as.character(lepi.annot.trees.family[[i]]$tip.label)
  #}

get.cherries <- function (phy) if (is.na(phy[1])) NA else sapply(which(node.depth(phy) == 2), function (x)  phy$tip.label[getDescendants(phy, x)], simplify = F)

lepi_tree_lineages_all <- bind_rows(lepi_tree_lineage_ts)

lepi_lineage_pairs_final <- list()
outgroup.trees <- list()

lepi_tribes_pairs <- sapply(lepi.annot.trees.tribes, get.cherries, simplify = F)
lepi_tribes_pairs_named <- sapply(seq_along(lepi_tribes_pairs), function (x) sapply(lepi_tribes_pairs[[x]], function (y) sapply(y, function (z) sapply(z, function (q) (lepi_tree_lineage_ts[[x]] %>% filter(Tip == q) %>% pull(Tribe))[1])), simplify=F), simplify=F)
for (tre in seq_along(lepi_tribes_pairs_named))
{

  for (pair in seq_along(lepi_tribes_pairs_named[[tre]]))
  {
    if (all(lepi_tree_lineage_ts[[tre]] %>% filter(Tribe %in% lepi_tribes_pairs_named[[tre]][[pair]]) %>% pull(TribeIsMono)) &
        !(length(lepi_tree_lineage_ts[[tre]] %>% filter(Tribe %in% lepi_tribes_pairs_named[[tre]][[pair]]) %>% pull(Subfamily) %>% unique)) > 1)
    {

      lepi_lineage_pairs_final <- c(lepi_lineage_pairs_final, list(lepi_tribes_pairs_named[[tre]][[pair]]))

      pairmrca <- getMRCA(lepi.annot.trees.tribes[[tre]], unlist(lepi_tribes_pairs[[tre]][[pair]]))
      parent.node <- getParent(lepi.annot.trees.tribes[[tre]], pairmrca)
      outgroup.parent.node <- lepi.annot.trees.tribes[[tre]]$edge[lepi.annot.trees.tribes[[tre]]$edge[,1] == parent.node, 2]
      
      outgroup.descs <- getDescendants(lepi.annot.trees.tribes[[tre]], 
                                       outgroup.parent.node[outgroup.parent.node != pairmrca]
      )

      poss.outgroups <- lepi.annot.trees.tribes[[tre]]$tip.label[
        outgroup.descs[
          outgroup.descs <= Ntip(lepi.annot.trees.tribes[[tre]])
          ]
        ]
      
      print(lepi.annot.trees.tribes[[tre]]$tip.label)
      poss.outgroups.tree <- keep.tip(lepi.annot.trees.tribes[[tre]], c(poss.outgroups, lepi_tribes_pairs[[tre]][[pair]]))
      
      poss.outgroups.tree$tip.label <- sapply(poss.outgroups.tree$tip.label, function (x) lepi_tree_lineage_ts[[tre]] %>% filter(Tip == x) %>% pull(Tribe))
      outgroup.trees <- c(outgroup.trees, list(poss.outgroups.tree))
    }
    
  }
}

lepi_subfam_pairs <- sapply(lepi.annot.trees.subfam, get.cherries, simplify = F)
lepi_subfam_pairs_named <- sapply(seq_along(lepi_subfam_pairs), function (x) sapply(lepi_subfam_pairs[[x]], function (y) sapply(y, function (z) sapply(z, function (q) (lepi_tree_lineage_ts[[x]] %>% filter(Tip == q) %>% pull(Subfamily))[1])), simplify=F), simplify=F)
for (tre in seq_along(lepi_subfam_pairs_named))
{
  for (pair in seq_along(lepi_subfam_pairs_named[[tre]]))
  {
    if (all(lepi_tree_lineage_ts[[tre]] %>% filter(Subfamily %in% lepi_subfam_pairs_named[[tre]][[pair]]) %>% pull(SubfamIsMono)) &
        !(length(lepi_tree_lineage_ts[[tre]] %>% filter(Subfamily %in% lepi_subfam_pairs_named[[tre]][[pair]]) %>% pull(Family) %>% unique)) > 1)
    {
      lepi_lineage_pairs_final <- c(lepi_lineage_pairs_final, list(lepi_subfam_pairs_named[[tre]][[pair]]))
      
      pairmrca <- getMRCA(lepi.annot.trees.subfam[[tre]], unlist(lepi_subfam_pairs[[tre]][[pair]]))
      parent.node <- getParent(lepi.annot.trees.subfam[[tre]], pairmrca)
      outgroup.parent.node <- lepi.annot.trees.subfam[[tre]]$edge[lepi.annot.trees.subfam[[tre]]$edge[,1] == parent.node, 2]
      
      outgroup.descs <- getDescendants(lepi.annot.trees.subfam[[tre]], 
                                       outgroup.parent.node[outgroup.parent.node != pairmrca]
      )
      
      poss.outgroups <- lepi.annot.trees.subfam[[tre]]$tip.label[
        outgroup.descs[
          outgroup.descs <= Ntip(lepi.annot.trees.subfam[[tre]])
          ]
        ]
      
      poss.outgroups.tree <- keep.tip(lepi.annot.trees.subfam[[tre]], c(poss.outgroups, lepi_subfam_pairs[[tre]][[pair]]))
      
      poss.outgroups.tree$tip.label <- sapply(poss.outgroups.tree$tip.label, function (x) lepi_tree_lineage_ts[[tre]] %>% filter(Tip == x) %>% pull(Subfamily))
      outgroup.trees <- c(outgroup.trees, list(poss.outgroups.tree))
    }
    
  }
}

lepi_family_pairs <- sapply(lepi.annot.trees.family, get.cherries, simplify = F)
lepi_family_pairs_named <- sapply(seq_along(lepi_family_pairs), function (x) sapply(lepi_family_pairs[[x]], function (y) sapply(y, function (z) sapply(z, function (q) (lepi_tree_lineage_ts[[x]] %>% filter(Tip == q) %>% pull(Family))[1])), simplify=F), simplify=F)
for (tre in seq_along(lepi_family_pairs_named))
{
  for (pair in seq_along(lepi_family_pairs_named[[tre]]))
  {
    if (all(lepi_tree_lineage_ts[[tre]] %>% filter(Family %in% lepi_family_pairs_named[[tre]][[pair]]) %>% pull(FamIsMono)) &
        !(length(lepi_tree_lineage_ts[[tre]] %>% filter(Family %in% lepi_family_pairs_named[[tre]][[pair]]) %>% pull(Superfamily) %>% unique)) > 1)
    {
      lepi_lineage_pairs_final <- c(lepi_lineage_pairs_final, list(lepi_family_pairs_named[[tre]][[pair]]))
      
      pairmrca <- getMRCA(lepi.annot.trees.family[[tre]], unlist(lepi_family_pairs[[tre]][[pair]]))
      parent.node <- getParent(lepi.annot.trees.family[[tre]], pairmrca)
      outgroup.parent.node <- lepi.annot.trees.family[[tre]]$edge[lepi.annot.trees.family[[tre]]$edge[,1] == parent.node, 2]
      
      outgroup.descs <- getDescendants(lepi.annot.trees.family[[tre]], 
                                       outgroup.parent.node[outgroup.parent.node != pairmrca]
      )
      
      poss.outgroups <- lepi.annot.trees.family[[tre]]$tip.label[
        outgroup.descs[
          outgroup.descs <= Ntip(lepi.annot.trees.family[[tre]])
          ]
        ]
      
      poss.outgroups.tree <- keep.tip(lepi.annot.trees.family[[tre]], c(poss.outgroups, lepi_family_pairs[[tre]][[pair]]))
      
      poss.outgroups.tree$tip.label <- sapply(poss.outgroups.tree$tip.label, function (x) lepi_tree_lineage_ts[[tre]] %>% filter(Tip == x) %>% pull(Family))
      outgroup.trees <- c(outgroup.trees, list(poss.outgroups.tree))
    }
    
  }
}

lepi.lepindex.raw <- read.csv("Data\\Lepi_MajorLineages\\Lepi_MajorLineages_Lepindex.csv", stringsAsFactors = F)

lepi.valid.spp <- lepi.lepindex.raw %>% distinct(Name_data, .keep_all = T)%>% filter(str_detect(Name_data, "Valid Name")) %>% group_by(Genus) %>% summarise(Valid_spp = n() - 1)

lepi.hosts <- read.csv("Data\\Lepi_MajorLineages\\Lepi_MajorLineages_HOSTS.csv", stringsAsFactors = F)
lepi.hosts$Genus <- toupper(lepi.hosts$Genus)
lepi.hosts.bysp <- lepi.hosts %>% group_by(Species) %>% summarise(Family=first(Family), Genus=first(Genus), Host_spp = length(unique(Host_name)))

angiosperm.tree <- read.nexus("Data\\Janssens_et_al_2020_Angiosperm_Tree.new")

lepi_lineage_pairs_m <- sapply(lepi_lineage_pairs_final, function (x) sapply(x, function (y) y))
rownames(lepi_lineage_pairs_m) <- NULL

lepi.lineages.spp <- apply(lepi_lineage_pairs_m, c(1,2), function (x) lepi.valid.spp %>% filter(Genus %in% (lepi_taxonomy %>% filter(x == Family | x == Subfamily | x == Tribe) %>% pull(Genus))) %>% pull(Valid_spp) %>% sum)

lepi.lineages.hosts <- apply(lepi_lineage_pairs_m, c(1,2), function (x) length(lepi.hosts %>% filter(Genus %in% (lepi_taxonomy %>% filter(x == Family | x == Subfamily | x == Tribe) %>% pull(Genus))) %>% pull(Host_family) %>% unique()))
lepi.lineages.hosts.spp <- apply(lepi_lineage_pairs_m, c(1,2), function (x) length(lepi.hosts %>% filter(Genus %in% (lepi_taxonomy %>% filter(x == Family | x == Subfamily | x == Tribe) %>% pull(Genus))) %>% pull(Host_name) %>% unique()))
lepi.lineages.hosts.pgen <- apply(lepi_lineage_pairs_m, c(1,2), function (x) sum((lepi.hosts.bysp %>% filter(Genus %in% (lepi_taxonomy %>% filter(x == Family | x == Subfamily | x == Tribe) %>% pull(Genus))) %>% pull(Host_spp)) > 1))/lepi.lineages.hosts.spp
lepi.lineages.hosts.mean <- apply(lepi_lineage_pairs_m, c(1,2), function (x) mean((lepi.hosts.bysp %>% filter(Genus %in% (lepi_taxonomy %>% filter(x == Family | x == Subfamily | x == Tribe) %>% pull(Genus))) %>% pull(Host_spp))))
lepi.lineages.hosts.pd <- apply(lepi_lineage_pairs_m, c(1,2), function (x) sum(keep.tip(angiosperm.tree, intersect(angiosperm.tree$tip.label, gsub(" ","_", lepi.hosts %>% filter(Genus %in% (lepi_taxonomy %>% filter(x == Family | x == Subfamily | x == Tribe) %>% pull(Genus))) %>% pull(Host_name) %>% unique)))$edge.length))

lepi.lineages.out <- rbind(lepi.lineages.spp, 
                            lepi.lineages.hosts, 
                            lepi.lineages.hosts.spp, 
                            lepi.lineages.hosts.pgen, 
                            lepi.lineages.hosts.mean, 
                            lepi.lineages.hosts.pd, 
                            log(lepi.lineages.spp), 
                            log(lepi.lineages.hosts), 
                            log(lepi.lineages.hosts.spp), 
                            asin(sqrt(lepi.lineages.hosts.pgen)), 
                            log(lepi.lineages.hosts.mean), 
                            log(lepi.lineages.hosts.pd)
                            )

lepi.lineages.out <- t(rbind(1:dim(lepi.lineages.out)[2], 
                              lepi_lineage_pairs_m, 
                              lepi.lineages.out, 
                              abs(lepi.lineages.out[13,] - lepi.lineages.out[14,]), 
                              (lepi.lineages.out[15,] - lepi.lineages.out[16,]) * sign(lepi.lineages.out[1,] - lepi.lineages.out[2,]),
                              (lepi.lineages.out[17,] - lepi.lineages.out[19,]) * sign(lepi.lineages.out[1,] - lepi.lineages.out[2,]), 
                              (lepi.lineages.out[19,] - lepi.lineages.out[20,]) * sign(lepi.lineages.out[1,] - lepi.lineages.out[2,]), 
                              (lepi.lineages.out[21,] - lepi.lineages.out[22,]) * sign(lepi.lineages.out[1,] - lepi.lineages.out[2,]), 
                              (lepi.lineages.out[23,] - lepi.lineages.out[24,]) * sign(lepi.lineages.out[1,] - lepi.lineages.out[2,])
                              )
                        )
colnames(lepi.lineages.out) <- c("Pair_no", 
                                  "Lineage_first", 
                                  "Lineage_last", 
                                  "N_spp_first", 
                                  "N_spp_last", 
                                  "Host_families_first", 
                                  "Host_families_last", 
                                  "Host_species_first", 
                                  "Host_species_last", 
                                  "Hosts_pgen_first", 
                                  "Hosts_pgen_last", 
                                  "Hosts_mean_first", 
                                  "Hosts_mean_last", 
                                  "Hosts_pd_first", 
                                  "Hosts_pd_last", 
                                  "log_N_spp_first", 
                                  "log_N_spp_last", 
                                  "log_Host_families_first", 
                                  "log_Host_families_last", 
                                  "log_Host_species_first", 
                                  "log_Host_species_last", 
                                  "asin_sqrt_Hosts_pgen_first", 
                                  "asin_sqrt_Hosts_pgen_last", 
                                  "log_Hosts_mean_first", 
                                  "log_Hosts_mean_last", 
                                  "log_Hosts_pd_first", 
                                  "log_Hosts_pd_last", 
                                  "N_spp_std", 
                                  "Host_families_std", 
                                  "Host_species_std", 
                                  "Hosts_pgen_std", 
                                  "Hosts_mean_std", 
                                  "Hosts_pd_std")
lepi.lineages.out.t <- as_tibble(lepi.lineages.out) %>% mutate(Rank = ifelse(str_detect(Lineage_first, "inae"), "Subfamily", (ifelse(str_detect(Lineage_first, "idae"), "Family", "Tribe")))) %>% filter(N_spp_first != 0, N_spp_last != 0, Host_families_first != 0, Host_families_last != 0) %>% distinct(Lineage_first, .keep_all = T)

write.csv(lepi.lineages.out.t, "Lepi_MajorLineages_Contrasts.csv")

lineages.hosts.fam.lm <- summary(lm(data=lepi.lineages.out.t, formula=as.numeric(N_spp_std) ~ as.numeric(Host_families_std)-1))$coeff
lineages.hosts.spp.lm <- summary(lm(data=lepi.lineages.out.t %>% filter(is.finite(as.numeric(Host_species_std))), formula=as.numeric(N_spp_std) ~ as.numeric(Host_species_std)-1))$coeff
lineages.hosts.pgen.lm <- summary(lm(data=lepi.lineages.out.t %>% filter(is.finite(as.numeric(Hosts_pgen_std))), formula=as.numeric(N_spp_std) ~ as.numeric(Hosts_pgen_std)-1))$coeff
lineages.hosts.mean.lm <- summary(lm(data=lepi.lineages.out.t %>% filter(is.finite(as.numeric(Hosts_mean_std))), formula=as.numeric(N_spp_std) ~ as.numeric(Hosts_mean_std)-1))$coeff
lineages.hosts.pd.lm <- summary(lm(data=lepi.lineages.out.t %>% filter(is.finite(as.numeric(Hosts_pd_std))), formula=as.numeric(N_spp_std) ~ as.numeric(Hosts_pd_std)-1))$coeff

lineages.hosts.fam.pl <- ggplot(lepi.lineages.out.t) + geom_point(aes(y = as.numeric(N_spp_std), x = as.numeric(Host_families_std))) + ylab("Lineage contrast in clade size") +  xlab("Lineage contrast in number of host plant families")+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5), text = element_text(family="serif")) + geom_abline(linetype = "dashed", size=1, colour="red", slope=lineages.hosts.fam.lm[1], intercept = 0)
lineages.hosts.spp.pl <- ggplot(lepi.lineages.out.t %>% filter(is.finite(as.numeric(Host_species_std)))) + geom_point(aes(y = as.numeric(N_spp_std), x = as.numeric(Host_species_std))) + ylab("Lineage contrast in clade size") +  xlab("Lineage contrast in number of host plant species")+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5), text = element_text(family="serif")) + geom_abline(linetype = "dashed", size=1, colour="grey80", slope=lineages.hosts.spp.lm[1], intercept = 0)
lineages.hosts.pgen.pl <- ggplot(lepi.lineages.out.t %>% filter(is.finite(as.numeric(Hosts_pgen_std)))) + geom_point(aes(y = as.numeric(N_spp_std), x = as.numeric(Hosts_pgen_std))) + ylab("Lineage contrast in clade size") +  xlab("Lineage contrast in proportion of host generalist species")+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5), text = element_text(family="serif")) + geom_abline(linetype = "dashed", size=1, colour="grey80", slope=lineages.hosts.pgen.lm[1], intercept = 0)
lineages.hosts.mean.pl <- ggplot(lepi.lineages.out.t %>% filter(is.finite(as.numeric(Hosts_mean_std)))) + geom_point(aes(y = as.numeric(N_spp_std), x = as.numeric(Hosts_mean_std))) + ylab("Lineage contrast in clade size") +  xlab("Lineage contrast in mean host species per Lepidopteran species")+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5), text = element_text(family="serif")) + geom_abline(linetype = "dashed", size=1, colour="red", slope=lineages.hosts.mean.lm[1], intercept = 0)
lineages.hosts.pd.pl <- ggplot(lepi.lineages.out.t %>% filter(is.finite(as.numeric(Hosts_pd_std)))) + geom_point(aes(y = as.numeric(N_spp_std), x = as.numeric(Hosts_pd_std))) + ylab("Lineage contrast in clade size") +  xlab("Lineage contrast in mean host species per Lepidopteran species")+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5), text = element_text(family="serif")) + geom_abline(linetype = "dashed", size=1, colour="red", slope=lineages.hosts.pd.lm[1], intercept = 0)

lineages.cordat <- cor(lepi.lineages.out.t %>% select_at(vars(contains("std"), -contains("rate"), -contains("time"))) %>% mutate_all(function (x) as.numeric(simplify2array(x))) %>%   filter(is.finite(Hosts_pd_std), !is.na(Host_families_std)))
lineages.cordat[upper.tri(lineages.cordat)] <- NA
ggplot(melt(lineages.cordat), aes(x=Var1, y=Var2)) + geom_tile(aes(fill=value, text=value)) + geom_text(aes(Var1, Var2, label = round(value,2)), size = 5) + scale_color_gradient2(na.value="grey80", low="red", mid="white", high="skyblue", midpoint=0, aesthetics="fill", name="Pearson", limit=c(-1,1))  + theme(axis.text.x = element_text(hjust=1, angle=45))

