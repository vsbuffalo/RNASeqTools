## AllGenerics.R

setGeneric("mdsPlot", signature="x", function(x, conds=NULL, cex=1, ...) standardGeneric("mdsPlot"))
setGeneric("joinColumns", signature="x", function(x, groups, fun=rowSums) standardGeneric("joinColumns"))
setGeneric("sumTrancripts2Genes", signature="x", function(x, mapping) standardGeneric("sumTrancripts2Genes"))
