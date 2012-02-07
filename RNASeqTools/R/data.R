## data.R -- functions for manipulating RNA-Seq data

setMethod("joinColumns", "data.frame",
          function(x, groups) {
            # the result will have the same number of rows as x, and the same
            # number of columns as unique integers in groups.
            if (ncol(x) != length(groups))
              stop("groups must be the same length as the number of columns of x.")
            out <- vector('list', length(unique(groups)))
            
            out <- lapply(unique(groups), function(g) {
              ci <- which(groups == g)
              rowSums(x[, ci])
            })
            names(out) <- unique(groups)
            d <- as.data.frame(do.call(cbind, out))
            
            rownames(d) <- rownames(x)
            d
          })

setMethod("joinRows", "data.frame", 
          function(x, mapping) {
            if (!(is.data.frame(x) && is.data.frame(mapping)))
              stop("x and mapping must be dataframes.")
            if (ncol(mapping) != 2)
              stop("The mapping dataframe must be two columns.")
            
            # reorder mapping file such that it corresponds to rownames of x (transcripts)
            ii <- match(rownames(x), mapping[, 2])
            not.found <- sum(is.na(ii))
            
            if (not.found > 0)
              warning(sprintf("%d row IDs were not found in the mapping file. If this number is high, check the the second column contains the correct IDs (often transcripts or exon IDs)!", not.found))
            
            tmp <- split(x, list(genes=mapping[ii, 1]))
            do.call(rbind, lapply(tmp, colSums))
          })

setMethod("joinRows", "matrix",
          function(x, mapping) {
            joinRows(as.data.frame(x), mapping)
          })

setMethod("joinColumns", "matrix",
          function(x, groups) {
            joinColumns(as.data.frame(x), groups)
          })

## setMethod("pairwiseComparison", "CountDataSet",
##           )



standardDESeq <- function(x, comparison) {
  cds <- newCountDataSet(x, comparison)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only")
  res <- nbinomTest(cds, comparison[1], comparison[2])
  res$comparison <- paste(comparison, collapse=";")
  res
}

setMethod("allComparisons", "factor",
          function(x, conds) {
            stopifnot(nlevels(conds) == 2)
            comparisons <- as.matrix(expand.grid(x[conds == levels(conds)[1]],
                                       x[conds == levels(conds)[2]]))
            comparisons
          })

setMethod("allComparisons", "character",
          function(x, conds) {
            tmp <- as.factor(conds)
            allComparisons(x, conds)
          })

runComparisons <- function(x, comparisons, analysis.function=standardDESeq,
                               db.name="toptable.db", mc.cores=2) {
  out <- mclapply(1:nrow(comparisons),
                  function(i) standardDESeq(x[, comparisons[i, ]], comparisons[i, ]), mc.cores=2)
  out
}
            
combineComparisons <- function(x) do.call(rbind, x)

queryComparisons <- function(x, p.cutoff=0.1, cutoff.col='padj') {
  interesting.i <- x[, cutoff.col] <= p.cutoff
  tmp <- sort(table(x$id[interesting.i]), decreasing=TRUE)
  list(n.interesting=table(x$comparison[interesting.i]), interesting.genes=tmp[tmp > 1])
}


setMethod("rowVars", "data.frame", function(x) apply(x, 1, var))
setMethod("rowVars", "matrix", function(x) apply(x, 1, var))

setMethod("rowCV", "data.frame", function(x) apply(x, 1, var)/rowMeans(x))
setMethod("rowCV", "matrix", function(x) apply(x, 1, var)/rowMeans(x))
setMethod("rowSCV", "data.frame", function(x) rowCV(x)^2)
setMethod("rowSCV", "matrix", function(x) rowCV(x)^2)

setMethod("varFilter", "data.frame", function(x, cutoff=0.4) {
  genevars <- rowVars(x)
  varcutoff <- quantile(genevars, probs=c(cutoff))
  i <- genevars >= varcutoff
  x[i, ]
})


