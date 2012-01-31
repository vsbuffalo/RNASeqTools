## AllGenerics.R

setGeneric("mdsPlot", signature="x", function(x, conds=NULL, cex=1, ...) standardGeneric("mdsPlot"))
setGeneric("joinColumns", signature="x", function(x, groups) standardGeneric("joinColumns"))
setGeneric("plotGeneDistribution", signature="x", function(x, showTopGenes=FALSE, showZeros=FALSE) standardGeneric("plotGeneDistribution"))
setGeneric("joinRows", signature="x", function(x, mapping) standardGeneric("joinRows"))
setGeneric("countDensity", signature="x", function(x) standardGeneric("countDensity"))
setGeneric("plotLibSizeSensitivity", signature="x", function(x) standardGeneric("plotLibSizeSensitivity"))
setGeneric("MAplot", signature="x", function(x, conds=conds, pval=pval, highlight=highlight,
                       smear=smear, interact=interact, xlog=xlog, xlab=xlab,
                       ylab=ylab, cex=cex, returnSmear=returnSmear, adj.fun=adj.fun)
           standardGeneric("MAplot"))
