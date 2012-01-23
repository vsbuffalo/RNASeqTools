\name{joinRows-methods}
\docType{methods}
\alias{joinRows}
\alias{joinRows-methods}
\alias{joinRows,data.frame-method}
\alias{joinRows,matrix-method}
\title{Take Transcript-level Counts and Sum to Gene-level Counts}
\description{
  \code{joinRows} takes a dataframe of gene expression data
  (rows are genes, columns are samples) and a dataframe of mapping
  information (i.e. transcript or exon  membership to genes) in which the
  first column contains gene IDs and the second column contains
  transcript/exon IDs.
}

\usage{
  joinRows(x, mapping)
}

\arguments{
  \item{x}{a dataframe of gene expression data.}
  \item{mapping}{a two column dataframe, in which the first column
  contains gene IDs and the second column contains transcript/exon IDs.}
}

\value{A dataframe in which the rows of transcripts from the same genes
  have been summed.}

\examples{
  library(pasilla)
  data(pasillaExons)

  ## create a mapping file 
  d <- counts(pasillaExons)
  tmp <- do.call(rbind, strsplit(rownames(counts(pasillaExons)), ":"))
  gene.exon.mapping <- data.frame(gene=tmp[, 1], exon=rownames(counts(pasillaExons)))

  ## join exons to genes
  d.joined <- joinRows(d, gene.exon.mapping)
  stopifnot(all(d.joined["FBgn0000256", ] == colSums(d[grep("FBgn0000256", rownames(d)), ])))
}

\keyword{methods}
\keyword{}
