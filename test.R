## test.R -- load package, test features

## Load package
## library(RNASeqTools)
library(DESeq)
sapply(list.files("RNASeqTools/R/", pattern="\\.R$", full.names=TRUE), source)

## Load pasillaGenes
library(pasilla)
data(pasillaGenes)
d <- as.data.frame(counts(pasillaGenes))

