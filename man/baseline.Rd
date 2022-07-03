\name{baseline}
\alias{baseline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
baseline
}
\description{
baseline generates a spline to fit the wandering baseline of the PPG time series and subtracts it from the time series.
}
\usage{
baseline(inx, iny, o, dat, sp, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inx}{
A vector of all inflection point x-coordinates in the PPG time series
}
 \item{iny}{
A vector of all inflection point y-coordinates in the PPG time series
}
 \item{o}{
A vector of all O points
}
 \item{dat}{
The original PPG time series (preprocessed)
}
 \item{sp}{
The original PPG time series, polynomial spline form
}
 \item{plot}{
Logical, plots baseline spline and the original PPG time series
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{baseCor }{time series with baseline wander removed}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Given baseline drift is indicative of vasomotion and respiratory modulation, it may be of interest to extract the spline.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
baseCor <- baseline(inx = inflexX, iny = inflexY, o = o_orig,
                    dat = undetrended, sp = splinePoly, plot=F)
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
}
