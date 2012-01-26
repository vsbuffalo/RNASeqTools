\name{plotLibSizeSensitivity-methods}
\docType{methods}
\alias{plotLibSizeSensitivity-methods}
\alias{plotLibSizeSensitivity,CountDataSet-method}
\alias{plotLibSizeSensitivity,data.frame-method}
\alias{plotLibSizeSensitivity,matrix-method}
\title{Plot Size Factors Against Number Zero Counts}
\description{
  Library sizes vary naturally in RNA-Seq experiments. Small libraries
  are essentially sampling a cell or tissue's to a lesser extent than
  larger libraries. This could lead to fewer genes detection (it has
  been show that deeper sequencing always leads to more genes being
  found).

  This plot looks at the number of genes with zero counts (genes that
  were not detected in a sample) with the library size. Using total
  counts is not a robust estimate of library size, as a few highly
  expressed genes could dominate total sample counts. Instead,
  plotLibSizeSensitivity uses the size factors estimated by
  \code{DESeq}'s \code{estimateSizeFactorsForMatrix}. See
  \code{\link{estimateSizeFactorsForMatrix}} for the procedure.
}

\usage{
  plotLibSizeSensitivity(x)
}

\arguments{
  \item{x}{either a dataframe or an object of class \code{CountDataSet} from the
    \code{DESeq} package.}
}

\section{Methods}{
\describe{

  \item{\code{signature(x = "data.frame")}}{
    This makes a plot showing the relationship between library size
    and the number of genes/transcripts with zero counts.
  }

  \item{\code{signature(x = "matrix")}}{
    Same as above.
  }

  \item{\code{signature(x = "CountDataSet")}}{
    This uses a similar method as for \code{x} with class
    \code{data.frame}, except it extracts the condition information from
    the phenotype data and colors by treatment group.
  }
}}

\examples{

  library(pasilla)
  data(pasillaGenes)

  # colored by treatment groups
  plotLibSizeSensitivity(pasillaGenes)

  # color by samples (because this is a dataframe)
  d <- as.data.frame(counts(pasillaGenes))

  plotLibSizeSensitivity(d)
}

\seealso{estimateSizeFactorsForMatrix}
\keyword{methods}
\keyword{diagnostics}
