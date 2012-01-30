## plots.R -- functions to create RNA-seq related plots

setMethod("mdsPlot", "data.frame", 
          function(x, conds=NULL, cex=1, ...) {
            d <- dist(t(x))
            
            mds.fit <- cmdscale(d, eig=TRUE, k=2)
            
            mds.d <- data.frame(x1=mds.fit$points[, 1],
                                x2=mds.fit$points[, 2],
                                labels=colnames(x))
            if (!is.null(conds))
              mds.d$treatment <- as.factor(conds)
            
            if (!is.null(conds)) {
              p <- xyplot(x2 ~ x1, group=treatment, data=mds.d, panel=function(x, y, ..., groups, subscripts) {
                              
                panel.text(x, y, mds.d$labels[subscripts], cex=cex, col=trellis.par.get()$superpose.line$col[groups])
              }, ...)
            } else {
              p <- xyplot(x2 ~ x1, data=mds.d, panel=function(x, y, ..., groups, subscripts) {
                panel.text(x, y, mds.d$labels[subscripts], cex=cex)
              }, ...)
            }
            print(p)

            invisible(mds.d)
          })

setMethod("mdsPlot", "CountDataSet",
          function(x, conds=NULL, cex=1, ...) {
            if (all(is.na(sizeFactors(x))))
              x <- estimateSizeFactors(x)
            tmp <- estimateDispersions(x, method="blind")
            vsd <- getVarianceStabilizedData(tmp)
            if (is.null(conds)) {
              # try to extract conditions from phenotype data
              conds <- pData(x)$condition
            }
            mdsPlot(vsd, conds, cex, ...)
          })

setMethod("mdsPlot", "matrix",
          function(x, conds=NULL, cex=1, ...) mdsPlot(as.data.frame(x), conds, cex, ...))


setMethod("countDensity", "data.frame", function(x) {
  tmp <- stack(x)
  densityplot(~ values, group=ind, data=tmp, scales=list(x = list(log=TRUE)),
              xlab="counts (log10 scale)", ylab="density", plot.points=FALSE)
})

setMethod("countDensity", "matrix", function(x) {
  tmp <- as.data.frame(x)
  countDensity(tmp)
})

setMethod("countDensity", "CountDataSet", function(x) {  
  countDensity(as.data.frame(counts(x)))
})


geneDistribution <- function(x) {
  # Order genes from high counts to low counts, and then return a
  # dataframe with percent of genes contributing and percent of total
  # counts.
  i <- quantile(x, probs=seq(1, 0, -0.01))
  total <- sum(x)
  y <- sapply(i, function(p) sum(x[x >= p]))
  num.genes <- sapply(i, function(p) sum(x >= p))
  data.frame(x=1:length(y), y=y, y.percent.total=y/total, num.genes)
}
  
setMethod("plotGeneDistribution", "data.frame",
          function(x, showTopGenes=FALSE, showZeros=FALSE) {
            tmp <- apply(x, 2, geneDistribution)
            d <- do.call(rbind, tmp)
            names <- local({
              tmp <- do.call(rbind, strsplit(gsub("(.*)\\.([0-9]+)%", "\\1;;;\\2", rownames(d)), ";;;"))
              vars <- data.frame(sample=tmp[, 1], percent=as.numeric(tmp[, 2]))
              vars
            })
            d <- as.data.frame(cbind(d, names))
            rownames(d) <- NULL
            p <- xyplot(y.percent.total ~ x, group=sample, data=d, type='l',
                        xlab="% of genes contributing (ordered high to low)",
                        ylab="% of total counts",
                        auto.key=list(columns=3, space="bottom"),
                        par.settings=simpleTheme(pch=20),
                        panel=function(x, y, groups, ...) {
                          if (showTopGenes) {
                            tmp <- d[d$percent == 99, ]
                            panel.points(tmp$x, tmp$y.percent,
                                         col=trellis.par.get("superpose.line")$col[tmp$sample])
                            ltext(x=tmp$x+8, y=tmp$y.percent, labels=format(tmp$y, big.mark=","),
                                  col=trellis.par.get("superpose.line")$col[tmp$sample], 
                                  pos=1, offset=1, cex=0.8)
                          }
                          if (showZeros) {
                            tmp <- aggregate(d$y, list(sample=d$sample), max)
                            tmp <- d[d$sample==tmp$sample & d$y == tmp$x, ]
                            tmp <- aggregate(tmp$percent, list(sample=tmp$sample), max)
                            tmp <- d[tmp$sample == d$sample & tmp$x == d$percent, ]
                            panel.abline(v=tmp$x, lty=2, lwd=0.4,
                                         col=trellis.par.get("superpose.line")$col[tmp$sample])
                            
                          }
                          panel.xyplot(x, y, groups, ...)
                        })
            print(p)
            invisible(d)
          })

setMethod("plotGeneDistribution", "matrix", function(x, showTopGenes=FALSE, showZeros=FALSE)
          plotGeneDistribution(as.data.frame(x), showTopGenes, showZeros))

setMethod("plotGeneDistribution", "CountDataSet", function(x, showTopGenes=FALSE, showZeros=FALSE) {
  d <- counts(x)
  plotGeneDistribution(d, showTopGenes, showZeros)
})


 
### TO ADD: function that takes vector of p-values, gene lengths and looks for relationship
mart <- useMart("ensembl")
ensembl <- useDataset("dmelanogaster_gene_ensembl", mart)

plotLengthPval <- function(lengths, pvals, highlight=0.1) {
  plot(log10(lengths), -log10(pvals), type='n',
       xlab="gene or transcript length (log10)", ylab="-log10 p-value")
  points(log10(lengths), -log10(pvals), col=ifelse(pvals <= highlight, "red", "black"),
         pch=19, cex=0.3)
  smallest.dig.gene <- min(na.exclude(lengths[pvals <= highlight]))
  abline(v=log10(smallest.dig.gene), lty=2, lwd=0.4)
  invisible(smallest.dig.gene)
}


plotLengthSensitivity <- function(lengths, pvals, length.range=c(0, 7000), bin.size=100) {
  tmp <- data.frame(lengths, pvals)
  d <- subset(tmp, lengths >= length.range[1] & lengths <= length.range)

  lcuts <- cut2(d$lengths, m=bin.size)
  
  dd <- aggregate(d$pvals, list(length=lcuts), function(x) sum(na.exclude(x < 0.05))/length(na.exclude(x)))
  plot(dd)
}

setMethod("plotLibSizeSensitivity", "data.frame",
          function(x) {
            locfunc <- median ## see getMethod("estimateSizeFactors", "CountDataSet")
            size.factors <- estimateSizeFactorsForMatrix(x, locfunc)
            num.zeros <- apply(x==0, 2, sum)
            d <- data.frame(sample=colnames(x), size.factors, num.zeros)
            p <- xyplot(num.zeros ~ size.factors, group=sample, data=d,
                        pch=19, auto.key=list(columns=3, space="bottom"),
                        par.settings=simpleTheme(pch=20),
                        xlab="size factors (estimated using procedure in DESeq)",
                        ylab="number of genes/transcripts with zero counts",
                        main=sprintf("Library Size and Number of Zero Counts\n(Spearman Correlation is %.2f)", cor(num.zeros, size.factors, method="spearman")))
            print(p)
            invisible(d)
          })

setMethod("plotLibSizeSensitivity", "matrix",
          function(x) plotLibSizeSensitivity(as.data.frame(x)))

setMethod("plotLibSizeSensitivity", "CountDataSet",
          function(x) {
            conds <- pData(x)$condition
            locfunc <- median ## see getMethod("estimateSizeFactors", "CountDataSet")
            cds <- estimateSizeFactors(x)
            size.factors <- sizeFactors(cds)
            num.zeros <- apply(counts(x)==0, 2, sum)
            d <- data.frame(sample=colnames(counts(x)), condition=conds, size.factors, num.zeros)
            p <- xyplot(num.zeros ~ size.factors, group=condition, data=d,
                        pch=19, #auto.key=list(columns=3, space="bottom"),
                        par.settings=simpleTheme(pch=20),
                        xlab="size factors (estimated using procedure in DESeq)",
                        ylab="number of genes/transcripts with zero counts",
                        main=sprintf("Library Size and Number of Zero Counts\n(Spearman Correlation is %.2f)", cor(num.zeros, size.factors, method="spearman")),
                        panel=function(x, y, groups, subscripts, ...) {
                          panel.text(x, y, d$sample[subscripts],
                                     col=trellis.par.get()$superpose.line$col[groups])
                          panel.xyplot(x, y, groups, ...)
                        })
            print(p)
            invisible(d)
          })


## TO ADD: MA-plot with interact=true option for gene selection.
calcMA <- function(x, conds, xlog) {
  a <- xlog(rowMeans(x))
  unique.conds <- unique(conds)
  stopifnot(length(unique.conds) == 2)
  mean.a <- rowMeans(x[, unique.conds[1] == conds])
  mean.b <- rowMeans(x[, unique.conds[2] == conds])
  
  m <- log2(mean.a/mean.b)
  d <- data.frame(m, a)
  d
}

MAplot <- function(x, conds=NULL, pval=NULL, highlight=0.1, interact=FALSE, xlog=log10,
                   xlab="mean counts (log10)", ylab="log fold change") {

  d <- calcMA(x, conds, xlog)
  # Remove NAs (due to 0 counts), but plot them separately
  x.adj <- x[which(!is.finite(d$m)), ]
  x.adj <- x.adj + 1
  d.adj <- calcMA(x.adj, conds, xlog)

  if (!is.null(pval)) {
    d$pval <- pval
    plot(d$a, d$m, type='n', xlab=xlab, ylab=ylab)
    points(d$a, d$m, pch=19, cex=0.1, col=ifelse(d$pval <= highlight, "red", "black"))
    points(runif(nrow(d.adj), -1.5, -0.5), d.adj$m, pch=19, cex=0.1, col=ifelse(d$pval <= highlight, "green", "purple"))
    
  } else {
    plot(d$a, d$m, type='n', xlab=xlab, ylab=ylab)
    points(d$a, d$m, pch=19, cex=0.1)
    points(runif(nrow(d.adj), -1.5, -0.5), d.adj$m, pch=19, cex=0.1, col="purple")
  }
  
}
