## AllGenerics.R

setGeneric("mdsPlot", signature="x", function(x, conds=NULL, cex=1, ...) standardGeneric("mdsPlot"))
setGeneric("joinColumns", signature="x", function(x, groups) standardGeneric("joinColumns"))
setGeneric("plotGeneDistribution", signature="x", function(x, showTopGenes=FALSE, showZeros=FALSE) standardGeneric("plotGeneDistribution"))
setGeneric("joinRows", signature="x", function(x, mapping) standardGeneric("joinRows"))
