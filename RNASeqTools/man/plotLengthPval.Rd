\name{plotLengthPval}
\alias{plotLengthPval}
\title{Plot Gene/Transcript Length vs p-values}
\description{
  Longer genes/transcripts are more likely to be detected in RNA-Seq
  experiments. \code{plotLengthPval} takes a look at how strong this
  effect can be on a particular data set. Note that this does not
  necessarily mean there is technical bias against finding shorter
  genes/transcripts that are differentially expressed; there could be
  biological reasons why shorter genes/transcripts are differentially
  expressed at lower rates.
}
\usage{
  plotLengthPval(lengths, pvals, highlight = 0.1)
}

\arguments{
  \item{lengths}{a numerical vector of lengths of genes.}
  \item{pvals}{a numerical vector of p-values in the same order as the
  lengths vector.}
  \item{highlight}{a numerical value between 0 and 1 indicating which
  genes to highlight (color as red).}
}

\details{
  \code{plotLengthPval} creates a plot with the -log10 p-values on the
  y-axis and the log10 gene/transcript lengths on the
  x-axis. Significant genes (with p-value less than the value specified
  in \code{highlight}) are colored red. A vertical line indicates the
  length of the gene/transcript with the smallest length that is
  differentially expressed.
}

\value{
  The length of the differentially expressed gene/transcript with the
  smallest length is returned silently.
}

\author{
  Vince Buffalo <vsbuffalo@ucdavis.edu>
}

\seealso{
  \code{\link{plotGeneDistribution}}
}

\examples{
  library(pasilla)
  library(biomaRt)
  library(DESeq)
  data(pasillaGenes)

  ## Caution: this could take a while!

  ## extact count data and make a plot with raw counts
  d <- counts(pasillaGenes)

  mart <- useMart("ensembl")
  ensembl <- useDataset("dmelanogaster_gene_ensembl", mart)
  length.df <- getBM(attributes=c("ensembl_gene_id", "start_position",
    "end_position"), filters="ensembl_gene_id",
    values=rownames(counts(pasillaGenes)), mart=ensembl)
  length.df$length <- abs(length.df$end_position -
    length.df$start_position)

  cds <- newCountDataSet(counts(pasillaGenes),
    pData(pasillaGenes)$condition)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  res <- nbinomTest(cds, "treated", "untreated")
  res$length <- length.df$length[match(res$id, length.df$ensembl_gene_id)]
  plotLengthPval(res$length, res$padj)
}

\keyword{diagnostics}
