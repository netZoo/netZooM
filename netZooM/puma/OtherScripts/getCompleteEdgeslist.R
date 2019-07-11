### get a complete edgelist for a prior
  # NOTE: this doesn't work with priors with weighted edges.
    # the script would need to be adapted to handle such a prior

library(igraph)
library(reshape2)

prior <- read.delim("prior.txt", header=F) # substitute this with the motif prior you want to make a complete edgelist for. this should be an input file in PANDA/PUMA format: TF GENE edge

prior <- prior[order(as.character(prior[,2])),] # order the motif prior alphabetically
prior <- prior[order(as.character(prior[,1])),] # i added "as.character" to be sure that genes are ordered in the same way as the expression data (which also needs to be ordered based on characters). sometimes it happens that factors in dataframes are ordered differently from character vectors and this can create problems

g <- graph.data.frame(prior[,1:2]) # make igraph object
gra <- get.adjacency(g,sparse=FALSE) # convert into adjacency matrix. this is a complete matrix of TF+gene vs TF+gene

gra <- gra[1:length(unique(prior[,1])), (length(unique(prior[,1]))+1):ncol(gra)] # remove TF-TF or gene-gene interactions from this complete matrix

edg <- melt(gra) # convert to a "PANDA/PUMA" prior object
edg <- edg[order(as.character(edg[,2])),] # to be sure that things are ordered alphabetically, let's order again
edg <- edg[order(as.character(edg[,1])),] # might want to check if this is really necessary if you want to reduce computing time
write.table(edg, "prior.txt", sep="\t", quote=F, row.names=F, col.names=F)
