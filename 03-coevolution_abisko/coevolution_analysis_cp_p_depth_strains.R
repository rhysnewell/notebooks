# coevolution based on depth strains 

library("ape")
library("phylogram")
library("dendextend")
library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ For the no indels trees ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#read trees in parenthetic format (Newick or New Hampshire format)
cp.1.tree = read.tree(file="06-parsnp_trees/c_palustre_no_indels.tree")
cp.1.tree$tip.label = sapply(str_split(cp.1.tree$tip.label, "38[_/.]"), "[", 2)
tp.1.tree = read.tree(file="06-parsnp_trees/t_primus_no_indels.tree")
tp.1.tree$tip.label = sapply(str_split(tp.1.tree$tip.label, "28[_/.]"), "[", 2)

#compute the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths 
cp.1.distances = cophenetic(cp.1.tree)
tp.1.distances = cophenetic(tp.1.tree)

#Need to figure out how to make the HP matrix (rectangular matrix with hosts as rows and parasites as columns)
#The matrix contains 1's when a host-parasite link has been observed between the host (row) and parasite (column), and 0's otherwise.
hp1 = read.csv("06-parsnp_trees/hp1.headers.txt", sep="\t", header = TRUE, row=1)
hp1 = as.matrix(hp1)

#Test host-parasite coevolution 
parafit(tp.1.distances, cp.1.distances, hp1, nperm = 1000, test.links = TRUE, correction = "cailliez")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ For the indel trees ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#read trees in parenthetic format (Newick or New Hampshire format)
cp.2.tree = read.tree(file="cp_indels.tree")
tp.2.tree = read.tree(file="tp_indels.tree")

#compute the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths 
cp.2.distances = cophenetic(cp.2.tree)
tp.2.distances = cophenetic(tp.2.tree)

#Need to figure out how to make the HP matrix (rectangular matrix with hosts as rows and parasites as columns)
#The matrix contains 1's when a host-parasite link has been observed between the host (row) and parasite (column), and 0's otherwise.
hp2 = read.csv("hp2.headers.txt", sep="\t", header = TRUE, row=1)
hp2 = as.matrix(hp2)

#Test host-parasite coevolution 
parafit(tp.2.distances, cp.2.distances, hp2, nperm = 1000, test.links = TRUE, correction = "cailliez")



# ~~~~~~~~~~~~~~

#makes the trees ultrametric 
cp.ultra = chronos(cp.1.tree)
tp.ultra = chronos(tp.1.tree)

#roots the tree 
cp.ultra$root.edge = 0
tp.ultra$root.edge = 0

#make binary 
cp.binary = multi2di(cp.ultra)
tp.binary = multi2di(tp.ultra)

#ploy the things
plot.phylo(cp.binary, direction="leftwards")
plot.phylo(tp.binary)

#need to relabel the tips the same things for this to work (I think)
tanglegram(cp.binary, tp.binary)
