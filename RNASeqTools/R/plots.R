## plots.R -- functions to create RNA-seq related plots

setMethod("mdsPlot", "data.frame", 
          function(x, conds=NULL, cex=1, ...) {
            d <- dist(x)
            
            mds.fit <- cmdscale(d, eig=TRUE, k=2)
            
            mds.d <- data.frame(x1=mds.fit$points[, 1],
                                x2=mds.fit$points[, 2],
                                labels=rownames(x))
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
            return(p)
          })
