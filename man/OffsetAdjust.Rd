\name{OffsetAdjust}
\alias{OffsetAdjust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OffsetAdjust
}
\description{
OffsetAdjust identifies the most appropriate offset value to be used in undetrending an ISO study ppg time
series (for inputting into the Undetrend function). It is used within the FindUndetrendingParams function
once a factor value has been determined.
}
\usage{
OffsetAdjust(ppg3, ppg, u = UnDetrend, factor_value, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppg3}{
the ppg time series adjusted for a given factor value
}
  \item{ppg}{
the original ppg time series
}
  \item{u}{
Undetrend function
}
  \item{factor_value}{
the appropriate factor value for restoring proper morphology to the ppg time series
}
  \item{plot}{
logical, if TRUE, plots the entire time series for each iteration of offset adjustment, thereby
visualising the adjustment as it occurs
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{offset_value }{determines the overall trend of the time series}
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
OffsetAdjust(ppg3, ppg, u = UnDetrend, factor_value, plot = F)
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
