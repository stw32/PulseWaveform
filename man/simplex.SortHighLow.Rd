\name{simplex.SortHighLow}
\alias{simplex.SortHighLow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simplex.SortHighLow
}
\description{
simplex.SortHighLow calculates the vertices (parameter sets) of a simplex with the highest, second highest and lowest ChiSq
values (expressed as row numbers of the simplex matrix).
}
\usage{
simplex.SortHighLow(vec_ChiSq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec_ChiSq}{
a vector of goodness of fit (ChiSq) values corresponding to each row (/parameter set) of the simplex
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{low,nHigh,high }{vector of the lowest, second highest, and highest ChiSq values in the simplex}
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
extrema <- simplex.SortHighLow(chiSq)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
