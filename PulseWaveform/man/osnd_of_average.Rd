\name{osnd_of_average}
\alias{osnd_of_average}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
osnd_of_average
}
\description{
osnd_of_average identifies OSND points on the average waveform and individual waves.
}
\usage{
osnd_of_average(aw, dp, diff, sr, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aw}{
average waveform, or an individual waveform
}
  \item{dp}{
diastolic peak (as established by the diast_pk function)
}
  \item{diff}{
difference in x-axis values of the average waveform due to removal of NA values (outputted from diast_pk as xShift)
}
  \item{sr}{
sampling rate
}
  \item{plot}{
logical - for plotting the inputted waveform with identified OSND points
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{osnd }{OSND values for the inputted waveform}
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
osnd <- osnd_of_average(avWave, dp = dPeak, diff = 0, sr = samplingRate, plot = F)
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
