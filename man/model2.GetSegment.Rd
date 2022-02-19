\name{model2.GetSegment}
\alias{model2.GetSegment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model2.GetSegment
}
\description{
model2.GetSegment extracts an original waveform from the preprocessed ppg time series. It is used by a number
of other functions.
}
\usage{
model2.GetSegment(ppg, limits)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppg}{
ppg time series
}
  \item{limits}{
limits of segment of time series to be extracted
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{result }{segment of ppg time series corresponding to a single waveform}
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
model2.GetSegment(ppg, limits)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
