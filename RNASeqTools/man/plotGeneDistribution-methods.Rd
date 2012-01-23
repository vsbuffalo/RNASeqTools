\name{plotGeneDistribution-methods}
\docType{methods}
\alias{plotGeneDistribution}
\alias{plotGeneDistribution-methods}
\alias{plotGeneDistribution,data.frame-method}
\alias{plotGeneDistribution,matrix-method}
\alias{plotGeneDistribution,CountDataSet-method}
\title{Plot Gene Distribution}
\description{
  \code{plotGeneDistribution} is a function to plot the percent of genes
  contribution to total lane counts, for each sample. This is a useful
  diagnostic to see if a few genes dominate total lane counts. Options
  are available to show the counts of the top 1\% of genes and the
  position of the genes with zero counts in the distribution.
}

\usage{
  plotGeneDistribution(x, showTopGenes=FALSE, showZeros=FALSE)
}

\arguments{
  \item{x}{either a dataframe or an object of class \code{CountDataSet} from the
    \code{DESeq} package.}
  \item{showTopGenes}{a logical indicating whether to show the
    position and raw counts of the top 1\% genes.}
  \item{showZeros}{a logical indicating whether to draw a vertical line
    at the position in which genes have zero counts.}
}



\section{Methods}{
  \describe{

    \item{\code{signature(x = "data.frame")}}{
      Plot a gene count distribution for each sample.
    }

    \item{\code{signature(x = "data.frame")}}{
      Same as with a \code{data.frame}, but with \code{CountDataSet} input.
    }

    \item{\code{signature(x = "matrix")}}{
      Same as with a \code{data.frame}, but with \code{matrix} input.
    }
}}


\value{
  A dataframe produced before the plot is made, containing columns for
  percent of genes contributing (\code{x}), counts of these genes
  (\code{y}), counts as percent sample total (\code{y.percent.total}), a
  factor vector of samples (\code{sample}), and percent (\code{percent}).
}

\examples{
  library(pasilla)
  data(pasillaGenes)

  ## extact count data and make a plot with raw counts
  d <- counts(pasillaGenes)
  plotGeneDistribution(d, showTopGenes=TRUE, showZeros=TRUE)

  ## or, use CountDataSet directly
  plotGeneDistribution(pasillaGenes)
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}

\keyword{methods}
\keyword{diagnostics}
