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
