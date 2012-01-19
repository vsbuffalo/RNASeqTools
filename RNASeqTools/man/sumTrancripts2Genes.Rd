\name{sumTrancripts2Genes}
\alias{sumTrancripts2Genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sumTrancripts2Genes(x, mapping)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mapping}{
%%     ~~Describe \code{mapping} here~~
}
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
function (x, mapping) 
{
    if (!(is.data.frame(x) && is.data.frame(mapping))) 
        stop("x and mapping must be dataframes.")
    if (ncol(mapping) != 2) 
        stop("The mapping dataframe must be two columns.")
    ii <- match(rownames(x), mapping[, 2])
    not.found <- sum(is.na(ii))
    if (not.found > 0) 
        warning(sprintf("\%d transcripts were not found in the mapping file. If this number is high, check the the second column is the transcript IDs!", 
            not.found))
    tmp <- split(x, list(genes = mapping[ii, 1]))
    do.call(rbind, lapply(tmp, colSums))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
