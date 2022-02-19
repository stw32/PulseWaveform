\name{model2.SubtractExcessPeak}
\alias{model2.SubtractExcessPeak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model2.SubtractExcessPeak
}
\description{
 model2.SubtractExcessPeak removes modelled component sine waves from excess segments of the ppg time series.
}
\usage{
model2.SubtractExcessPeak(time, residue, peakParams)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
x-values of ppg time series
}
  \item{residue}{
excess component of ppg time series (before or after the extraction of other component waves)
}
  \item{peakParams}{
vector of three parameters to construct the component sine wave to be removed
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{result }{excess time series with modelled sine wave subtracted (residue)}
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
residue <- sep(data[,1],residue,par[8:10])
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
