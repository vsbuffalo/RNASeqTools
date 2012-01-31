\name{MAplot}
\alias{MAplot}
\alias{MAplot-methods}
\alias{MAplot,data.frame-method}
\alias{MAplot,CountDataSet-method}
\title{Plot Average Log Counts over Log Fold Change}
\description{\code{MAplot} produces a plot in which the x-axis contains
  the average counts across treatments (log-transformed) and the y-axis
  contains the log fold change between the treatment groups
  (log2-transformed). This type of plot was originally used for
  microarray data; this particular function is only suited for RNA-seq
  count data.

  This function contains many options to help the biologist or
  bioinformatician explore gene expression data. If one treatment group
  has all zero counts, its average will be zero giving a non-finite or
  not a number (NaN) log fold change. Normally, R would not plot these
  values. \code{MAplot} handles this by transforming them with the
  including function \code{addOne} specified to the \code{adj.fun}
  argument, which simply adds 1 to both  counts. This small amount of
  bias allows for easier interpretation. Other adjustment functions can
  be specified.

  The adjusted values can either be plot in the original MA plot, but
  the default behavior is to create a "smear" plot, with the adjusted
  value's log fold changes being plotted away from the normal value's MA
  plot.

  \code{MAplot} also allows user interaction through the
  \code{interact=TRUE} arguement.
}
\usage{
MAplot(x, conds = NULL, pval = NULL, highlight = 0.1, smear = TRUE, interact = FALSE, xlog = log10, xlab = "mean counts (log10)", ylab = "log fold change", cex = 0.2, returnSmear = FALSE, adj.fun = addOne)
}

\arguments{
  \item{x}{an object of class \code{data.frame} or \code{CountDataSet}
  that contains RNA-seq expression data.}
  \item{conds}{a factor, character, or integer vector indicating
  treatment membership of the columns of \code{x}.}
  \item{pval}{a numeric vector of p-values that correspond to the data
  in \code{x} (same order; this is the default in \code{DESeq}).}
  \item{highlight}{a numeric value. P-values less than or equal to this
  value will be colored red in the MA-plot.}
  \item{smear}{a logical indicating whether to "smear" non-finite
  log-fold changes outside of the MA-plot, or whether to plot the
  adjusted values and just color them differently.}
  \item{interact}{a logical indicating whether to allow the user to
  select points and print their labels (i.e. gene/transcript name). If
  \code{TRUE} the return output will be a character vector of the
  selected points.}
  \item{xlog}{a function in the log family that transforms the x-axis values.}
  \item{xlab}{a character value to be used as the x-axis label.}
  \item{ylab}{a character value to be used as the y-axis label.}
  \item{cex}{a numeric value that scales the point size.}
  \item{returnSmear}{a logical value. When \code{MAplot} smears points,
  it chooses their position at random. When \code{interact=FALSE},
  \code{MAplot} returns a dataframe of the plotted points; with
  \code{returnSmear=FALSE}, \code{MAplot} returns the original points
  (which have non-finite log fold changes). The option to return the
  random points with \code{returnSmear=TRUE} allows users to access the
  data as it was plotted.}
  \item{adj.fun}{a function that tells \code{MAplot} how to adjust
  counts such that they have finite log values. By default this is
  \code{addOne}, which simply adds 1 to counts that would lead to
  non-finite log fold changes.}
}

\section{Methods}{
\describe{

  \item{\code{signature(x = "data.frame")}}{
    This is most common method signature; a dataframe of counts (these
    should be normalized in most cases!) is passed through to
    \code{MAplot}.
  }

  \item{\code{signature(x = "matrix")}}{
    Same as above.
  }

  \item{\code{signature(x = "CountDataSet")}}{
    This uses a similar method as for \code{x} with class
    \code{data.frame}, except it extracts the \emph{normalized} counts
    from the \code{CountDataSet} object.
  }
}}


\value{
  If \code{interact=FALSE}, \code{MAplot} returns a dataframe of all the
  plotted points, with the very important except that the adjusted
  points in the smear which are given a random uniform position on the
  x-axis will be returned with their original (NaN or non-finite)
  values, unless \code{returnSmear=TRUE}. The \code{returnSmear} option was
  added to allow users to retrieve the exact points of the plot if
  needed. The random points are stored within the function, as they are
  needed for when \code{interact=TRUE}.

  If \code{interact=TRUE}, \code{MAplot} will return a vector of the
  gene/transcript names selected by the user, in the order in which they
  were selected.  
}

\author{
Vince Buffalo <vsbuffalo@ucdavis.edu>
}

\examples{
  library(pasilla)
  data(pasillaGenes)

  conds <- pData(pasillaGenes)$condition
  cds <- newCountDataSet(counts(pasillaGenes), conds)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersion(cds)
  res <- nbinomTest(cds, "treated", "untreated")

  MAplot(counts(cds, normalized=FALSE), conds, pval=res$padj)
  MAplot(counts(cds, normalized=TRUE), conds, pval=res$padj)
  MAplot(counts(cds, normalized=TRUE), conds, smear=FALSE,
    pval=res$padj)

  ## Not plotted as it would have to wait for user interaction
  ## MAplot(counts(cds, normalized=TRUE), conds, pval=res$padj,
  ## interact=TRUE)
  
}


\keyword{diagnostics}

