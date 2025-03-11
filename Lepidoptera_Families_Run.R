library(phytools)
library(tidyverse)
library(matrixcalc)
library(reshape2)

source("R\\baseml.extract.tree.R")
source("R\\my.drop.tip.R")
source("R\\ordered.ls.dir.R")
source("R\\phylo.average.brlen.R")
source("R\\welch.test.R")

lepi.fam.aln <- as.character(read.FASTA("Data\\Lepi_Families\\SA3_nucleotide_supermatrix_before_degen.fas"))
lepi.fam.parts <- read.table("Data\\Lepi_Families\\SA3_nucleotide_gene_boundaries.txt")
lepi.fam.parts.t <- lepi.fam.parts %>% mutate(first = sapply(V3, function (x) as.integer(str_split(x, '-')[[1]][1])), second = sapply(V3, function (x) as.integer(str_split(gsub(';', '', x), '-')[[1]][2]))) %>% select(-c(V2, V3))
lepi.fam.genes <- mapply(function (x,y) sapply(lepi.fam.aln, function (z) z[x:y], simplify = F), lepi.fam.parts.t[,2], lepi.fam.parts.t[,3], SIMPLIFY = F)
lepi.fam.covs <- sapply(lepi.fam.genes, function (x) sapply(x, function (y) sum(y %in% c('a','c','t','g'))))
lepi.fam.seqnames <- read.csv("Data\\Lepi_Families\\Lepi_Families_seqnames.csv")
colnames(lepi.fam.covs) <- sapply(str_split(lepi.fam.parts.t[,1], "\\."), function (x) x[[1]])
rownames(lepi.fam.covs) <- lepi.fam.seqnames$seqname
names(lepi.fam.aln)<- lepi.fam.seqnames$seqname
lepi.fam.ngenes <- (lepi.fam.covs > 0) %*% rep(1, dim(lepi.fam.covs)[2])

lepi.fam.S12.tre <- read.tree("Data\\Lepi_Families\\Fig_S12_Newick.tre")
get.cherries <- function (phy) sapply(which(node.depth(phy) == 2), function (x) phy$tip.label[getDescendants(phy, x)], simplify = F)

lepi.fam.rogue.tips <- c("Bombycoidea_Sphingidae_Hemileucinae_Hemileuca_maia_BuckH",
                        "Pyraloidea_Pyralidae_Spilomelinae_Spilomela_sp_INSerlTAIRAAPEI__75",
                        "Sesioidea_Castniidae_Castniinae_Synemon_plana_Cul",
                        "Sesioidea_Castniidae_Castniinae_Telchin_licus_SRR1204999",
                        "Palaephatoidea_Palaephatidae_Ptyssoptera_sp_Ptys",
                        "Palaephatoidea_Palaephatidae_Metaphatus_ochraceus_Met4",
                        "Palaephatoidea_Palaephatidae_Palaephatus_nielseni_Pro43D",
                        "Palaephatoidea_Palaephatidae_Palaephatus_luteolus_Pal")
lepi.fam.norogue.tre <- drop.tip(lepi.fam.S12.tre, lepi.fam.rogue.tips)
lepi.fam.families <- unique(sapply(str_split(lepi.fam.norogue.tre$tip.label, "_"), function (x) x[2]))
lepi.fam.family.tips <- sapply(lepi.fam.families, function (x) lepi.fam.norogue.tre$tip.label[str_detect(lepi.fam.norogue.tre$tip.label, x)], simplify = F)
lepi.fam.family.tree <- keep.tip(lepi.fam.norogue.tre, sapply(lepi.fam.family.tips, function (x) x[1]))
lepi.fam.family.tree$tip.label <- names(lepi.fam.family.tips)
lepi.fam.family.pairs <- get.cherries(lepi.fam.family.tree)

lepi.fam.pairs.sampled <- list()

for (pair in lepi.fam.family.pairs)
{
  sislengths  <- sapply(pair, function (x) length(lepi.fam.family.tips[[x]]))
  bigsister <- lepi.fam.family.tips[[pair[which.max(sislengths)]]]
  littlesister <- lepi.fam.family.tips[[pair[-which.max(sislengths)]]]

  bigsister.ngenes <- lepi.fam.ngenes[bigsister,]
  bigsister.sampled <- bigsister[order(bigsister.ngenes, decreasing = T)][1:min(sislengths)]

  #bigsister.sampled <- sample(bigsister, min(sislengths))

  lepi.fam.pairs.sampled <- c(lepi.fam.pairs.sampled,
                             list(list(bigsister.sampled, littlesister)))
}

outgroup.ids <- character(length(lepi.fam.pairs.sampled))

for (i in 1:length(lepi.fam.pairs.sampled))
{
  pairmrca <- getMRCA(lepi.fam.norogue.tre, unlist(lepi.fam.pairs.sampled[[i]]))
  parent.node <- getParent(lepi.fam.norogue.tre, pairmrca)
  outgroup.parent.node <- lepi.fam.norogue.tre$edge[lepi.fam.norogue.tre$edge[,1] == parent.node, 2]
  outgroup.descs <- getDescendants(lepi.fam.norogue.tre,
                                   outgroup.parent.node[outgroup.parent.node != pairmrca]
                                   )
  poss.outgroups <- lepi.fam.norogue.tre$tip.label[
    outgroup.descs[
      outgroup.descs <= Ntip(lepi.fam.norogue.tre)
      ]
    ]

  outgroup.ids[i] <- poss.outgroups[which.max(lepi.fam.ngenes[poss.outgroups,])]

}
 
# lepi.fam.pair.alns <- mapply(function (x,y) lepi.fam.aln[c(unlist(x),y)],
#                             lepi.fam.pairs.sampled,
#                             outgroup.ids,
#                             SIMPLIFY = F)
 


# lepi.fam.pair.alns.ungap <- list()
# lepi.fam.parts.ungap <- list()
# lepi.fam.pair.trees <- sapply(lepi.fam.pair.alns,
#                              function (x) keep.tip(lepi.fam.norogue.tre, names(x)),
#                                                    simplify =F)
 
# complete.codons <- function (x) sapply(x,
#                                        function (y) sapply(seq(1, length(y), 3),
#                                                            function (z) all(y[z:(z+2)] %in% c('a','t','c','g','A','C','T','G')
#                                                                             )
#                                                            )
#                                        )
# 
# for (aln in lepi.fam.pair.alns)
# {
#   aln.m <- sapply(aln, complete.codons)
#   n <- (length(aln) - 1)/2
#   print(paste0("Trimming alignment of ", dim(aln.m)[2] - 1, " sequences plus outgroup."))
#   aln.goodcols <- as.logical((((((aln.m %*%
#                          direct.sum(direct.sum(rep(1,n), rep(1,n)), 1)) > 0) %*%
#                        rep(1,3)
#                      )== 3) %x% rep(1,3))[,1])
# 
#   output.aln <- sapply(aln, function (x) x[aln.goodcols])
# 
#   lepi.fam.pair.alns.ungap <- c(lepi.fam.pair.alns.ungap, list(output.aln))
# 
#   aln.parts <- mapply(function (x,y) sum(aln.goodcols[x:y]),
#                              lepi.fam.parts.t$first,
#                              lepi.fam.parts.t$second)
# 
#   lepi.fam.parts.ungap <- c(lepi.fam.parts.ungap,
#                            list(tibble(locus=lepi.fam.parts.t$V1, bp=aln.parts) %>%
#                                 filter(bp != 0) %>%
#                                 mutate(terminus = cumsum(bp)))
#   )
# 
# }

# for (i in seq_along(lepi.fam.pair.alns.ungap)) 
# {
#   write.FASTA(as.DNAbin(t(lepi.fam.pair.alns.ungap[[i]])), 
#               file = paste0("Pair_alignments\\Lepi_Families\\lepi.fam_", 
#                             i, 
#                             ".fasta")
#               )
#   print(paste0("Writing alignment ", i, ", ", 
#         dim(lepi.fam.pair.alns.ungap[[i]])[1], " bp, ", 
#         dim(lepi.fam.pair.alns.ungap[[i]])[2], " seq"))
#   
#   write.csv(lepi.fam.parts.ungap[[i]], 
#             paste0("Data\\Lepi_Families\\lepi_fam_parts_", i, ".csv")
#             )
#   
#   write.tree(lepi.fam.pair.trees[[i]], 
#              file=paste0("Pair_trees\\Lepi_Families\\lepi_fam_topology_", i, ".tre")
#   )
# }

lepi.fam.baseml.trees <- sapply(ordered.ls.dir("PAML_outputs\\Lepi_Families", regexp = ".*baseml.txt"), baseml.extract.tree, simplify = F)
lepi.fam.code.trees <- sapply(ordered.ls.dir("PAML_outputs\\Lepi_Families", regexp = ".*codeml.txt"), codeml.extract.trees, simplify = F)

lepi_taxonomy <- read.csv("Data\\Lepi_Families\\Lepi_Families_taxonomy.csv", stringsAsFactors = F)

lepi.fam.lepindex.raw <- read.csv("Data\\Lepi_Families\\Lepi_Families_Lepindex_Raw.csv", stringsAsFactors=F)

trichopt.valid.spp <- read.csv("Data\\Lepi_Families\\Trichoptera_Global_Checklist.csv") %>% group_by(Family) %>% summarise(Valid_spp = sum(Valid_spp))
lepi.fam.valid.spp <- lepi.fam.lepindex.raw %>% distinct(Name_data, .keep_all = T)%>% filter(str_detect(Name_data, "Valid Name")) %>% group_by(Family) %>% summarise(Valid_spp = n() - 1)
lepi.fam.valid.spp <- bind_rows(lepi.fam.valid.spp, trichopt.valid.spp, data.frame(Family="Dryadaulidae", Valid_spp=50)) %>% arrange(Family)

lepi.fam.blens <- mapply(function (x,y,q) sapply(c(list(q),y[2:3]), function (z) sapply(x, function (a) phylo.average.brlen(my.drop.tip(z, z$tip.label[str_detect(paste(a, collapse=' '), z$tip.label, negate = T)], root.edge = 1)))), lepi.fam.pairs.sampled, lepi.fam.code.trees, lepi.fam.baseml.trees)
lepi.fam.famnames <- sapply(lepi.fam.pairs.sampled, function (x) c(str_split(x[[1]], "_")[[1]][2], str_split(x[[2]], "_")[[1]][2]))
trichopt.spp.contrasts <- read.csv("Data\\Lepi_Families\\Trichoptera_Global_Checklist.csv") %>% group_by(Family) %>% summarise(Valid_spp = sum(Valid_spp))
lepi.fam.spp.contrasts <- apply(lepi.fam.famnames, c(1,2), function (x) lepi.fam.valid.spp %>% filter(Family %in% x) %>% pull(Valid_spp))
lepi.fam.out <- t(rbind(1:26, c(rep("Trichoptera", 3), rep("Lepidoptera", 23)), lepi.fam.famnames, lepi.fam.blens, lepi.fam.spp.contrasts))
colnames(lepi.fam.out) <- c("Pair_no", "Order", "Family_first", "Family_last", "blen_first", "blen_last", "dN_first", "dN_last", "dS_first", "dS_last", "N_spp_first", "N_spp_last")
lepi.fam.out.t <- as_tibble(lepi.fam.out) %>% mutate_at(vars(blen_first, blen_last, dN_first, dN_last, dS_first, dS_last, N_spp_first, N_spp_last), as.numeric, .keep="unused") %>% mutate_all(simplify2array)

lepi.fam.pair.ages <- sapply(lepi.fam.pairs.sampled, 
                          function (x) node.depth.edgelength(
                            keep.tip(lepi.fam.norogue.tre, 
                                     unlist(x))
                          )[1]
)

angiosperm.tree <- read.nexus("Data\\Janssens_et_al_2020_Angiosperm_Tree.new")

lepi.fam_hosts <- read.csv("Data\\Lepi_Families\\Lepi_Families_HOSTS.csv", stringsAsFactors = F)
lepi.fam_hosts_final <- lepi.fam_hosts %>% filter(Family !='') %>% mutate(Host_family = sapply(str_split(Host_family, ' '), function (x) x[1]), Host_genus = sapply(str_split(Host_name, " "), function (x) toupper(x[1])))

lepi.fam_hosts_summ1 <- lepi.fam_hosts_final %>% group_by(Family,Species) %>% summarise(N_host_sp = length(unique(Host_name))) %>% summarise(Hosts_pgen = sum(N_host_sp > 1, na.rm=T)/n(), Hosts_mean = mean(N_host_sp, na.rm=T))
lepi.fam_hosts_summ2 <- lepi.fam_hosts_final %>% group_by(Family) %>% summarise(Host_families = length(unique(Host_family)), Host_species = length(unique(Host_name)))
lepi.fam_hosts_summ <- lepi.fam_hosts_summ1 %>% left_join(lepi.fam_hosts_summ2, by = "Family")

lepi.fam_hosts_notthere <- unlist(c(lepi.fam.out.t$Family_first, lepi.fam.out.t$Family_last))[!(unlist(c(lepi.fam.out.t$Family_first, lepi.fam.out.t$Family_last)) %in% lepi.fam_hosts_summ$Family)]
lepi.fam_hosts_notthere.t <- tibble(Family=lepi.fam_hosts_notthere, Hosts_pgen=NA, Hosts_mean=NA, Host_families=NA, Host_species=NA)
lepi.fam_hosts_summ_final <- bind_rows(lepi.fam_hosts_summ, lepi.fam_hosts_notthere.t)

lepi.fam_hosttrees_pd <- apply(lepi.fam.famnames, c(1,2), function (x) keep.tip(angiosperm.tree, intersect(angiosperm.tree$tip.label, gsub(" ","_", lepi.fam_hosts_final %>% filter(toupper(Genus) %in% (lepi_taxonomy %>% filter(x == Family) %>% pull(Genus))) %>% pull(Host_name) %>% unique))))

lepi.fam_famnames_pd <- apply(lepi.fam_hosttrees_pd, c(1,2), function (x) if(is.null(x[[1]])) NA else phylo.average.brlen(x[[1]]))

lepi.fam_hosts_out <- t(rbind(apply(lepi.fam.famnames, c(1,2), function (x) lepi.fam_hosts_summ_final %>% filter(Family == x) %>% pull(Host_families)), apply(lepi.fam.famnames, c(1,2), function (x) lepi.fam_hosts_summ_final %>% filter(Family == x) %>% pull(Host_species)), apply(lepi.fam.famnames, c(1,2), function (x) lepi.fam_hosts_summ_final %>% filter(Family == x) %>% pull(Hosts_pgen)), apply(lepi.fam.famnames, c(1,2), function (x) lepi.fam_hosts_summ_final %>% filter(Family == x) %>% pull(Hosts_mean)), lepi.fam_famnames_pd))
colnames(lepi.fam_hosts_out) <- c("Host_families_first", "Host_families_last", "Host_species_first", "Host_species_last", "Hosts_pgen_first", "Hosts_pgen_last", "Hosts_mean_first", "Hosts_mean_last", "Hosts_pd_first", "Hosts_pd_last")
lepi.fam_hosts_out.t <- as_tibble(lepi.fam_hosts_out)
lepi.fam.out.t <- bind_cols(lepi.fam.out.t, lepi.fam_hosts_out.t) %>% transmute_all(simplify2array)

lepi.fam.out.t$Age=lepi.fam.pair.ages       

lepi.fam.hosts.out.unsat.t <- lepi.fam.out.t %>% filter(dS_first < 2, dS_last < 2, Order != "Trichoptera")
lepi.fam.hosts.out.final.t <- lepi.fam.hosts.out.unsat.t %>% filter(!is.na(N_spp_first), !is.na(N_spp_last)) 

write.csv("Lepi_Families_Contrasts.csv")