\name{find_o}
\alias{find_o}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
find_o
}
\description{
find_o identifies the origin of the systolic peak for each beat. In the absence of a clear inflection point at the origin, O points are derived from inflection points in the first derivative.
}
\usage{
find_o(wx, inx, iny, d1p, sp)
}
\arguments{
  \item{wx}{
A vector of x-coordinates corresponding to w points on the PPG time series
}
  \item{inx}{
A vector of all inflection point x-coordinates in the PPG time series
}
  \item{iny}{
A vector of all inflection point y-coordinates in the PPG time series
}
  \item{d1p}{
The first derivative time series, polynomial spline form
}
  \item{sp}{
The original PPG time series, polynomial spline form
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{o }{vector of o points}
\item{inx }{vector of updated inflection point x-coordinates}
\item{iny }{vector of updated inflection point y-coordinates}
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
tmp <- find_o(wx = w$wX, inx = inflexX, iny = inflexY, d1p = deriv1Poly, sp = splinePoly)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
