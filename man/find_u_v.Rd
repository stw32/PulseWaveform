\name{find_u_v}
\alias{find_u_v}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
find_u_v
}
\description{
find_u_v identifies the points on the systolic upstroke that correspond to the two half maximum values on the corresponding first  derivative peak. The first, before w, is denoted U. The second, after w, is denoted V.
}
\usage{
find_u_v(wx, wy, d1, d1p, spline, sr, plot)
}
\arguments{
  \item{wx}{
A vector of x-coordinates corresponding to w points on the PPG time series
}
  \item{wy}{
A vector of y-coordinates corresponding to w points on the PPG time series
}
  \item{d1}{
The first derivative time series, discrete form
}
  \item{d1p}{
The first derivative time series, polynomial spline form
}
  \item{spline}{
The original PPG time series, polynomial spline form
}
  \item{sr}{
samplingRate
}
  \item{plot}{
Logical, will plot the first derivative time series with U and V values as points if set to true
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{uX}{A vector of x-xoordinates corresponding to U points on the PPG time series}
\item{uY}{A vector of y-xoordinates corresponding to U points      "       "        }
\item{vX}{A vector of x-xoordinates corresponding to V points      "       "        }
\item{vY}{A vector of y-xoordinates corresponding to V points      "       "        }
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
find_u_v(wx = w$wX, wy = w$wY, d1 = deriv1, d1p = deriv1Poly,
     spline = splinePoly, sr = samplingRate, plot=F)
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
