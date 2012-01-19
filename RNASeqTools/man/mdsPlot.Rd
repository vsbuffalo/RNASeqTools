\name{mdsPlot}
\alias{mdsPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
\code{mdsPlot} takes a set of gene expression data in a dataframe and
plots a diagnostic MDS plot of the samples.
}
\usage{
mdsPlot(x, conds = NULL, cex = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a dataframe of gene expression results, with rows as genes and
    columns as samples.}
  \item{conds}{an optional character vector indicating the treatment groups of
    each sample.}
  \item{cex}{a numerical value giving the amount by which plotting text
    and symbols should be magnified relative to the default (see
    \code{par}) for more details.}
  \item{\dots}{other arguments passed to lattice's \code{xyplot}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, conds = NULL, cex = 1, ...) 
{
    d <- dist(x)
    mds.fit <- cmdscale(d, eig = TRUE, k = 2)
    mds.d <- data.frame(x1 = mds.fit$points[, 1], x2 = mds.fit$points[, 
        2], labels = rownames(x))
    if (!is.null(conds)) 
        mds.d$treatment <- as.factor(conds)
    if (!is.null(conds)) {
        p <- xyplot(x2 ~ x1, group = treatment, data = mds.d, 
            panel = function(x, y, ..., groups, subscripts) {
                panel.text(x, y, mds.d$labels[subscripts], cex = cex, 
                  col = trellis.par.get()$superpose.line$col[groups])
            }, ...)
    }
    else {
        p <- xyplot(x2 ~ x1, data = mds.d, panel = function(x, 
            y, ..., groups, subscripts) {
            panel.text(x, y, mds.d$labels[subscripts], cex = cex)
        }, ...)
    }
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
