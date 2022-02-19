\name{clean_wuv}
\alias{clean_wuv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
clean_wuv
}
\description{
clean_wuv first identifies erroneous beats and removes them, by identifying abnormal values of U and V, and the y-axis difference between them. It then calculates O-O intervals and inter-beat intervals, and removes beats preceding large intervals.

}
\usage{
clean_wuv(wuv, sp, inx, o, samp, bc, q = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wuv}{
a dataframe of all w, u and v values corresponding to detected and pre-cleaned beats in the time series
}
  \item{sp}{
the baseline-corrected PPG time series, polynomial spline form
}
  \item{inx}{
a vector of all inflection point x-coordinates in the PPG time series
}
  \item{o}{
a vector of all o points
}
  \item{samp}{
sampling rate
}
  \item{bc}{
the baseline-corrected PPG time series, discrete form
}
  \item{q}{
logical, will pause function and give the option of plotting rejected beats
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{ibi}{a vector of all inter-beat intervals (x-axis differences between sucessive W points)}
\item{oDiff}{a vector of all O-O intervals (x-axis differences between successive O points)}
\item{d }{a dataframe combining the following three structures:}
\item{wuv }{the inputted dataframe of w, u and v values, with rows corresponding to erroneous beats removed}
\item{diffVU}{a vector of y-axis differences between u and v for each beat - for scaling purposes}
\item{o2}{a vector of O points, with those corresponding to erroneous beats removed}
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
tmp <- clean_wuv(wuv = wuv, sp = splinePolyBC, inx = inflexX, o = o,
                 samp = samplingRate, bc = baseCor, q = F)
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
