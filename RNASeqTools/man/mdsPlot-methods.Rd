\name{mdsPlot-methods}
\docType{methods}
\alias{mdsPlot-methods}
\alias{mdsPlot,data.frame-method}
\title{Make a Diagnostic Plot of Samples from Gene Expression Data}
\description{
  \code{mdsPlot} takes a dataframe or an object of class
  \code{CountDataSet} from the \code{DESeq} and plots the samples
  according to multidimensional scaling.
}

\usage{
  mdsPlot(x, conds=NULL, cex=1, ...)
}

\arguments{
  \item{x}{either a dataframe or an object of class \code{CountDataSet} from the
    \code{DESeq} package.}
  \item{conds}{a character vector indicating the treatment groups of the
    samples (columns) in \code{x}.}
  \item{cex}{a numerical value giving the amount by which plotting text
    and symbols should be magnified relative to the default. See ?par
    for more details.}
  \item{\dots}{Further arguments passsed to lattice's \code{xyplot}.}
}
\keyword{methods}
\keyword{diagnostics,plots}
