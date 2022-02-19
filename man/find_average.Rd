\name{find_average}
\alias{find_average}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
find_average
}
\description{
find_average generates an average (mean) wave from a sample of segmented waveforms. Due to variance in waveform length, the tail section of the average wave is less robustly calculated, and relies on a degree of inference as to the likely trajectories of waveforms that terminate before the average length is reached.
}
\usage{
find_average(p, ao)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
dataframe consisting of all segmented waveforms
}
  \item{ao}{
a vector indicating the end point of each waveform
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{avWave }{a vector consisting of values for the average waveform}
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
average_wave <- find_average(p = pulse, ao = afterO)
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
