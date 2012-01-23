\name{countDensity-methods}
\docType{methods}
\alias{countDensity}
\alias{countDensity-methods}
\alias{countDensity,data.frame-method}
\alias{countDensity,matrix-method}
\alias{countDensity,CountDataSet-method}
\title{Plot Count Density}
\description{
  \code{countDensity} is a function to plot the density of counts for
  each sample, either from a dataframe or matrix, or a \code{CountDataSet}.
}

\usage{
  countDensity(x)
}

\arguments{
  \item{x}{either a dataframe or an object of class \code{CountDataSet} from the
    \code{DESeq} package.}
}


\section{Methods}{
  \describe{

    \item{\code{signature(x = "data.frame")}}{
      Plot density for each sample.
    }

    \item{\code{signature(x = "data.frame")}}{
      Same as with a \code{data.frame}, but with \code{CountDataSet} input.
    }

    \item{\code{signature(x = "matrix")}}{
      Same as with a \code{data.frame}, but with \code{matrix} input.
    }
}}


\examples{
  library(pasilla)
  data(pasillaGenes)

  ## extact count data and make a plot with raw counts
  d <- counts(pasillaGenes)
  countDensity(d)

  ## or, use CountDataSet directly
  countDensity(pasillaGenes)
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}

\keyword{methods}
\keyword{diagnostics}
