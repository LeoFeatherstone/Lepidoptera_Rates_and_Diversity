library(phytools)
library(tidyverse)
library(matrixcalc)
library(stringr)
library(fs)
library(reshape2)

source("R\\baseml.extract.tree.R")
source("R\\my.drop.tip.R")
source("R\\ordered.ls.dir.R")
source("R\\phylo.average.brlen.R")
source("R\\welch.test.R")

papi.aln <- as.character(read.dna("Data\\Papi_Genera\\S_2c_Dataset_with_COI.phy", format = "sequential"))
rownames(papi.aln) <- str_trim(rownames(papi.aln))
papi.parts <- read.csv("Data\\Papi_Genera\\Papi_Genera_partitions.csv")
papi.gaps <- matrix(papi.aln %in% c('a','c','g','t'), nrow = dim(papi.aln)[1])
rownames(papi.gaps) = rownames(papi.aln)
papi.covs <- mapply(function (x,y) papi.gaps[,x:y] %*% rep(1,y-x+1), papi.parts$first, papi.parts$last)
colnames(papi.covs) <- papi.parts$gid
rownames(papi.covs) <- rownames(papi.gaps)

papi.tre <- read.nexus("Data\\Papi_Genera\\S_3b_core_analysis_median_ages.tre")
papi.covs.complete <- papi.covs[(papi.covs[,-c(1,3,4,9)] > 0) %*% rep(1, dim(papi.covs)[2] - 4) == (dim(papi.covs)[2] - 4),-c(1,3,4,9)]
papi.covs.complete <- papi.covs.complete[-which(rownames(papi.covs.complete) == "RWH_96_0878_Macrosoma_bahiata"),]
papi.tre.complete <- drop.tip(papi.tre, setdiff(papi.tre$tip.label, str_trim(rownames(papi.covs.complete))))

get.cherries <- function (phy) sapply(which(node.depth(phy) == 2), function (x) phy$tip.label[getDescendants(phy, x)], simplify = F)

#papi.tip.pairs <- get.cherries(papi.tre.complete)

papi.tip.pairs.inc <- get.cherries(papi.tre)
papi.tip.pairs <- papi.tip.pairs.inc[sapply(papi.tip.pairs.inc, function (x) all(x %in% str_trim(rownames(papi.covs.complete))))]

papi.counts.complete <- papi.covs.complete %*% rep(1, dim(papi.covs.complete)[2])
names(papi.counts.complete) <- str_trim(names(papi.counts.complete))

outgroup.ids <- character(length(papi.tip.pairs))

for (i in 1:length(papi.tip.pairs))
{
  pairmrca <- getMRCA(papi.tre.complete, unlist(papi.tip.pairs[[i]]))
  parent.node <- getParent(papi.tre.complete, pairmrca)
  outgroup.parent.node <- papi.tre.complete$edge[papi.tre.complete$edge[,1] == parent.node, 2]
  outgroup.descs <- getDescendants(papi.tre.complete, 
                                   outgroup.parent.node[outgroup.parent.node != pairmrca]
  )
  poss.outgroups <- papi.tre.complete$tip.label[
    outgroup.descs[
      outgroup.descs <= Ntip(papi.tre.complete)
      ]
    ]
  outgroup.ids[i] <- names(which.max(papi.counts.complete[poss.outgroups,1]))
  
}

papi.pairs.with.outgroups <- mapply(function (x, y) c(x, y), papi.tip.pairs, outgroup.ids, SIMPLIFY=F)

papi.pair.ages <- sapply(papi.tip.pairs, function (x) node.depth.edgelength(keep.tip(papi.tre.complete, x))[1])
papi.pair.trees <- sapply(papi.pairs.with.outgroups, function (x) keep.tip(papi.tre.complete, x), simplify=F)

papi.good.poses <- unlist(apply(papi.parts[-c(1,3,4,9),], 1, function (x) x[2]:x[3]))
papi.pair.alns <- sapply(papi.pairs.with.outgroups, function (x) papi.aln[x, unlist(papi.good.poses)], simplify = F)


papi.pair.alns.ungap <- list()
  
for (aln in papi.pair.alns) {
  
  aln.m <- matrix(as.numeric(aln %in% c('a', 't', 'c', 'g')), nrow = dim(aln)[1] * 3)
  good.codons <- rep(1, dim(aln)[1] * 3) %*% aln.m
  good.sites <- good.codons %x% rep(1,3) == dim(aln)[1] * 3
  papi.pair.alns.ungap <- c(papi.pair.alns.ungap, list(aln[,good.sites]))
   
}

for (i in 1:length(papi.pair.alns.ungap)) {
  print(paste0("writing: Pair_alignments\\Papi_Genera\\papi_", 
                         i, 
                         ".fasta"))
  write.FASTA(as.DNAbin(papi.pair.alns.ungap[[i]]), 
              paste0("Pair_alignments\\Papi_Genera\\papi_", 
                     i, 
                     ".fasta")
              )

  write.tree(papi.pair.trees[[i]], 
             paste0("Pair_trees\\Papi_Genera\\papi_",
                    i,
                    ".tre")
             )
}

## SOURCE TO HERE - THEN USE TREES AND ALIGNMENTS ABOVE TO RUN BASEML & CODEML

papi.genera <- sapply(str_split(unlist(papi.tip.pairs), "_"), function (x) x[length(x[x!=""])-1])
papi.genus.names <- matrix(papi.genera, nrow=2)
rownames(papi.genus.names) <- c("Genus_first", "Genus_last")

papi.baseml.trees <- sapply(ordered.ls.dir("PAML_outputs\\Papi_Genera", regexp = ".*baseml.txt"), baseml.extract.tree, simplify = F)
papi.codeml.trees <- sapply(ordered.ls.dir("PAML_outputs\\Papi_Genera", regexp = ".*codeml.txt"), codeml.extract.trees, simplify = F)

papi.blens <- mapply(function (x,y) y$edge.length[y$edge[,2] <= Ntip(y)][y$tip.label %in% x][order(match(x, y$tip.label))], papi.tip.pairs, papi.baseml.trees)
papi.code.blens <- mapply(function (x,y) sapply(y[2:3], function (z) z$edge.length[z$edge[,2] <= Ntip(z)][z$tip.label %in% x][order(match(x, z$tip.label))]), papi.tip.pairs, papi.codeml.trees)

butmoth.lepindex.raw <- read.csv("Data\\Papi_Genera\\Lepi_Families_Lepindex_Raw.csv", stringsAsFactors=F)
butmoth.valid.spp <- butmoth.lepindex.raw %>% distinct(Name_data, .keep_all = T)%>% filter(str_detect(Name_data, "Valid Name")) %>% group_by(Genus) %>% summarise(Valid_spp = n() - 1)

papi.lepindex.raw <- read.csv("Data\\Papi_Genera\\Papi_Genera_Lepindex_Raw.csv", stringsAsFactors = F)
papi.valid.spp <- papi.lepindex.raw %>% distinct(Name_data, .keep_all = T)%>% filter(str_detect(Name_data, "Valid Name")) %>% group_by(Genus) %>% summarise(Valid_spp = n() - 1)
papi.valid.spp.final <- bind_rows(butmoth.valid.spp, papi.valid.spp) %>% distinct(Genus, .keep_all = T) %>% arrange(Genus)
papi.genera.notthere <- data.frame (Genus = toupper(papi.genera[!(toupper(papi.genera)  %in% papi.valid.spp.final$Genus)]), Valid_spp = NA)
papi.valid.spp.final <- bind_rows(papi.valid.spp.final, papi.genera.notthere)
papi.spp.contrasts <- apply(papi.genus.names, c(1,2), function (x) papi.valid.spp.final %>% filter(Genus == toupper(x)) %>% pull(Valid_spp))

angiosperm.tree <- read.nexus("Data\\Janssens_et_al_2020_Angiosperm_Tree.new")

papi.out <- t(rbind(matrix(1:length(papi.tip.pairs), nrow=1), papi.blens, papi.code.blens, papi.spp.contrasts, matrix(papi.pair.ages, nrow=1)))
colnames(papi.out) <- c("Pair_no", "blen_first", "blen_last", "dN_first", "dN_last", "dS_first", "dS_last", "N_spp_first", "N_spp_last", "Pair_age")
papi.out.t <- bind_cols(as_tibble(t(papi.genus.names)), as_tibble(papi.out))

papi_hosts <- read.csv("Data\\Papi_Genera\\Papi_Genera_HOSTS.csv", stringsAsFactors = F)
papi_hosts_final <- papi_hosts %>% filter(Family !='') %>% mutate(Host_family = sapply(str_split(Host_family, ' '), function (x) x[1]), Host_genus = sapply(str_split(Host_name, " "), function (x) toupper(x[1])))

papi_hosts_summ1 <- papi_hosts_final %>% group_by(Genus,Species) %>% summarise(Family=first(Family), N_host_sp = length(unique(Host_name))) %>% summarise(Family=first(Family), Hosts_pgen = sum(N_host_sp > 1, na.rm=T)/n(), Hosts_mean = mean(N_host_sp, na.rm=T))
papi_hosts_summ2 <- papi_hosts_final %>% group_by(Genus) %>% summarise(Host_families = length(unique(Host_family)), Host_species = length(unique(Host_name)))
papi_hosts_summ <- papi_hosts_summ1 %>% left_join(papi_hosts_summ2, by = "Genus")

papi_hosts_notthere <- unlist(papi.genus.names)[!(unlist(papi.genus.names) %in% papi_hosts_summ$Genus)]
papi_hosts_notthere.t <- tibble(Genus=papi_hosts_notthere, Family=NA, Hosts_pgen=NA, Hosts_mean=NA, Host_families=NA, Hosts_species=NA)
papi_hosts_summ_final <- bind_rows(papi_hosts_summ, papi_hosts_notthere.t)

papi_hosts_pd <- apply(papi.genus.names, c(1,2), function (x) sum(keep.tip(angiosperm.tree, intersect(angiosperm.tree$tip.label, gsub(" ","_", papi_hosts_final %>% filter(Genus == x) %>% pull(Host_name) %>% unique)))$edge.length))

papi.host.contrasts <- rbind(apply(papi.genus.names, c(1,2), function (x) papi_hosts_summ_final %>% filter(Genus == x) %>% pull(Hosts_pgen)), apply(papi.genus.names, c(1,2), function (x) papi_hosts_summ_final %>% filter(Genus == x) %>% pull(Hosts_mean)), apply(papi.genus.names, c(1,2), function (x) papi_hosts_summ_final %>% filter(Genus == x) %>% pull(Host_families)), apply(papi.genus.names, c(1,2), function (x) papi_hosts_summ_final %>% filter(Genus == x) %>% pull(Host_species)), papi_hosts_pd)
rownames(papi.host.contrasts)<- c("Hosts_pgen_first", "Hosts_pgen_last", "Hosts_mean_first", "Hosts_mean_last", "Host_families_first", "Host_families_last", "Host_species_first", "Host_species_last", "Hosts_pd_first", "Hosts_pd_last")

papi.hosts.out.t <- bind_cols(as_tibble(t(papi.genus.names)), as_tibble(papi.out.t), as_tibble(t(papi.host.contrasts)))
papi.hosts.out.final.t <- papi.hosts.out.t %>% filter(!is.na(N_spp_first), !is.na(N_spp_last))



### FINAL OUTPUT --------------------------------------------------------------------------------

write.csv(papi.hosts.out.final.t, "Papi_Genera_Contrasts.csv")