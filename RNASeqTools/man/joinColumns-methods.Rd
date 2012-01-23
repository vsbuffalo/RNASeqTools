\name{joinColumns-methods}
\docType{methods}
\alias{joinColumns}
\alias{joinColumns-methods}
\alias{joinColumns,data.frame-method}
\title{Join Columns of a Dataframe Together}
\description{
  \code{joinColumns} takes a dataframe and a vector of groups the same
  length of as the number of columns and joins
  the columns by adding them.
}

\usage{
  joinColumns(x, groups)
}

\arguments{
  \item{x}{a dataframe of gene expression data.}
  \item{groups}{a character or integer vector of groups with length
    equal to the number of columns in \code{x}. Elements with the same
  value will be joined by \code{fun}.}
}


\section{Methods}{
  \describe{
    \item{\code{signature(x = "data.frame")}}{
      This method joins the columns of dataframe \code{x} according to the
      groupings in \code{group}.
    }
  }
}

\value{
  A dataframe with the same number of columns as unique elements in
  \code{groups} and same number of rows as \code{x}.
}

\examples{
  ## Create some fake data, i.e. from single and paired-end mapping data
  ## that need to be summed.
  d <- data.frame(a_paired=rnorm(100), a_single=rnorm(100),
b_paired=rnorm(100), b_single=rnorm(100))

  ## This shows a good way to extract groups from colnames by removing
  ## single and paired strings from column names.
  groups <- gsub("([ab])_.*", "\\\\1", colnames(d)) ## alternate way

  d.joined <- joinColumns(d, groups)

  ## check that results are consistent.
  eps <- 1e4
  stopifnot(rowSums(d[, c(1, 2)]) - d.joined[, 1] < eps &&
    rowSums(d[, c(3, 4)]) - d.joined[, 2] < eps)
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\keyword{methods}
\keyword{data-manipulation}
