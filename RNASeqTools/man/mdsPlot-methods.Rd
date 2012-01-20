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

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\section{Methods}{
  \item{\code{signature(x = "data.frame")}}{
    This method makes a plot of the dimensionality-reduced data (using
    MDS) from a data.frame. It returns a dataframe of the data used to
    produced the plot.
  \item{\code{signature(x = "matrix")}}{
    Same as above.
  \item{\code{signature(x = "CountDataSet")}}{
    This method makes a plot of the dimensionality-reduced data (using
    MDS) from an object of class \code{CountDataSet}, the variance
    stabilized data is used.
}


\examples{
  library(pasilla)
  data(pasillaGenes)
  d <- counts(pasillaGenes)[, pData(pasillaGenes)$type=="paired-end"]
  conds <- as.character(pData(pasillaGenes)$condition[pData(pasillaGenes)$type=="paired-end"])
  mdsPlot(d, conds=conds)
}

\keyword{methods}
\keyword{diagnostics,plots}
