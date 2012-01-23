\name{mdsPlot-methods}
\docType{methods}
\alias{mdsPlot}
\alias{mdsPlot-methods}
\alias{mdsPlot,data.frame-method}
\alias{mdsPlot,matrix-method}
\alias{mdsPlot,CountDataSet-method}
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

\section{Methods}{
  \describe{
    \item{\code{signature(x = "data.frame")}}{
      This method makes a plot of the dimensionality-reduced data (using
      MDS) from a data.frame. It returns a dataframe of the data used to
      produced the plot.}
    \item{\code{signature(x = "matrix")}}{
      Same as above.}
    \item{\code{signature(x = "CountDataSet")}}{
      This method makes a plot of the dimensionality-reduced data (using
      MDS) from an object of class \code{CountDataSet}, the variance
      stabilized data is used.}
  }
}

\value{
  A dataframe containing the coordinates of the samples after the
  MDS procedure.
}

\examples{
  library(pasilla)
  data(pasillaGenes)

  ## extact count data and make a plot with raw counts
  d <- counts(pasillaGenes)
  conds <- as.character(pData(pasillaGenes)$condition)
  mdsPlot(d, conds=conds)

  ## same, but log10 transform data first and add 1 to avoid taking
  ## log10(0)
  mdsPlot(log10(d + 1), conds=conds)

  ## or, interact the object of class CountDataSet
  ## directly. This method first estimates the size
  ## factor, estimate the disperions using the "blind"
  ## method, and finally get the variance stabilized
  ## data and make an MDS plot. This method will
  ## extract conditions from the CountDataSet object.
  mdsPlot(pasillaGenes)
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}

\keyword{methods}
\keyword{diagnostics,plots}
