# File: baseml.extract.tree.R

# Quick functions to extract the output trees from baseml and codeml results via regular expression.
# Does no checking, so paml output files must currently be checked prior to use to ensure they have not terminated early or failed to run.

# Usage: 
# source("baseml.extract.tree.R")
# baseml_trees <- sapply(dir_ls(glob="*_baseml_output.txt"), baseml.extract.tree, simplify = F)
# class(baseml_trees[[1]]) // "phylo"
#
# codeml_trees <- sapply(dir_ls(glob="*_codeml_output.txt"), codeml.extract.trees, simplify = F)
# class (codeml_trees[[1]]$dS) //"phylo"
# class (codeml_trees[[1]]$dN) //"phylo"

require (ape)
require (fs)

# Extract trees from baseml output. Return is a list of phylo objects.
baseml.extract.tree <- function (path)
{

    paml <- readLines(path)
    rex <- regmatches(paml, regexec("^(.*);$", paml))
    tr <- read.tree(text = rev(rex[sapply(rex, function (x) length(x) > 0)])[[1]][1])

    return (tr)

}

# Extract trees from codeml output. Return is a list of lists, each with named members 
# $w (branch lengths are value of omega parameter for that branch), 
# $dS (branch lengths are synonymous substitutions/site), 
# $dN (branch lengths are nonsynonymous substitutions/site)
codeml.extract.trees <- function (path)
{

    paml <- readLines(path)
    rex <- regmatches(paml, regexec("^(.*);$", paml))
    trs <- sapply(1:3, function (x) read.tree(text = rev(
                                                     rex[sapply(rex, function (x) length(x) > 0)]
                                                 )[[x]][1]), simplify = F)
    names(trs) <- c("w", "dN", "dS")
    return (trs)

}
