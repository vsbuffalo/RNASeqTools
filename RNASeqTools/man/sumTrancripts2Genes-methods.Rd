\name{sumTrancripts2Genes-methods}
\docType{methods}
\alias{sumTrancripts2Genes-methods}
\alias{sumTrancripts2Genes,data.frame-method}
\title{Take Transcript-level Counts and Sum to Gene-level Counts}
\description{
  \code{sumTrancripts2Genes} takes a dataframe of gene expression data
  (rows are genes, columns are samples) and a dataframe of mapping
  information (mapping transcript membership to genes) in which the
  first column contains gene IDs and the second column contains
  transcript IDs.
}

\usage{
  sumTrancripts2Genes(x, mapping)
}

\arguments{
  \item{x}{a dataframe of gene expression data.}
  \item{mapping}{a two column dataframe, in which the first column
  contains gene IDs and the second column contains transcript IDs.}
}

\keyword{methods}
\keyword{}
