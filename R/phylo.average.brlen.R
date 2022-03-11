require (ape)

#tipmatched_mito_pairdiffs <- left_join(tipmatched_mito_data %>% group_by(quartet, pair_no, taxon, Genus, Species) %>% summarise_at(c("Length", "DepthRangeDeep.x", "lat_centroid"), mean) %>% summarise_at(c("Length", "DepthRangeDeep.x", "lat_centroid"), mean) %>% summarise_at(c("Length", "DepthRangeDeep.x", "lat_centroid"), mean) %>% arrange(pair_no), tipmatched_mito_data %>% select(taxon, valid_species), by = "taxon") %>% summarise_at(c("Length", "DepthRangeDeep.x", "lat_centroid", "valid_species"), function (x) abs(first(x) - last(x)))


#tipmatched_nuc_base <- ordered.ls.dir("/Users/aritchie/Dropbox/BromhamPostdoc/BigFish/4_rates/baseml", regexp="tipmatched_rag1.*txt")

#tipmatched_nuc_table <- read.csv("/Users/aritchie/Dropbox/BromhamPostdoc/BigFish/1_data/fishgenes_tipmatched_nuc_sampled.csv")

#tipmatched_nuc_base_tr <- sapply(tipmatched_nuc_base, baseml.extract.tree, simplify = F)

#tipmatched_nuc_table <- tipmatched_nuc_table %>% arrange(pair_no)

#trts <-  unlist(sapply(unique(tipmatched_nuc_table %>% pull(quartet)), function (i) sapply(unique(tipmatched_nuc_table %>% filter(quartet == i) %>% pull(pair_no)), function (x) phylo.average.brlen(keep.tip(tipmatched_nuc_base_tr[[i]], as.character(tipmatched_nuc_table %>% filter(pair_no == x) %>% pull(RAG1.acc)))), simplify= T)))

#trtb <- tibble(pair_no = tipmatched_nuc_table %>% pull(pair_no), av_blen <- trts)




phylo.average.brlen <- function (phy)
{

    po <- reorder.phylo(phy, order = "postorder")

    br <- array(0, dim = Nnode(phy) + Ntip(phy))

    for (i in 1:Nedge(phy))
    {
        br[po$edge[i, 1]] <- br[po$edge[i, 1]] + (po$edge.length[i] + br[po$edge[i, 2]])/2
    }

    root.edge.length = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
    singleton = ifelse(Ntip(phy) == 1, 2, 1) 

    return (br[Ntip(phy) + 1] * singleton + root.edge.length)

}

