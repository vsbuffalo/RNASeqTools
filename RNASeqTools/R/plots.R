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
  data.frame(x=1:length(y), y=y, y.percent.total=y/total)
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

