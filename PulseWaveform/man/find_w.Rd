\name{find_w}
\alias{find_w}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
find_w
}
\description{
find_w identifies peaks in the first derivative of the PPG time series (denoted w on the original PPG pulse wave).
A rolling window relative to heart rate is applied to identify beats and artefacts. Peaks identified within a window are confirmed as genuine / artefactual with a series of checks against thresholds derived from beats local to the window and across the entire time series.
}
\usage{
find_w(d1p, deriv1, sp, sr)
}

\arguments{
  \item{d1p}{
polynomial spline of first derivative of the PPG time series
}
  \item{deriv1}{
first derivative in discrete form
}
  \item{sp}{
polynomial spline of original time series
}
  \item{sr}{
sample rate
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{wX}{x axis coordinates of all w points}
\item{wY}{y axis coordinates of all w points on original PPG trace}
\item{wYD1}{y axis coordinates of all w points on 1st derivative trace}
}
\references{
Elgendi et al, 2018: Toward generating more diagnostic features from photoplethysmogram waveforms.
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
w <- find_w(d1p = deriv1Poly, deriv1 = deriv1, sp = splinePoly, sr = samplingRate)
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
