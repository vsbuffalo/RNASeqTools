## pairwise_test.R -- compare standard DESeq approach with pairwise
## approach
##
## Vince Buffalo <vsbuffalo@gmail.com>

library(DESeq)
library(RNASeqTools)
library(pasilla)
library(multicore)

data(pasillaGenes)

# original analysis
d <- counts(pasillaGenes)
conds <- pData(pasillaGenes)$condition
cds <- newCountDataSet(d, conds)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
res <- nbinomTest(cds, "treated", "untreated")
res <- res[order(res$padj), ]

# pairwise
pairs <- allComparisons(colnames(d), conds)
tmp <- combineComparisons(runComparisons(d, pairs))
out <- queryComparisons(tmp, p.cutoff=0.2)

# comparison of two approaches
p <- 1:200
n.intersect <- sapply(p, function(i) length(intersect(names(out$interesting.genes[1:i]), res$id[1:i])))
plot(p, n.intersect, type='l', col="green", 
     ylab="number of differentially expressed genes\nfound in both pairwise comparisons and replicated DESeq",
     xlab="top n differentially expressed genes by p-value")
abline(b=1, a=0, col="blue")
