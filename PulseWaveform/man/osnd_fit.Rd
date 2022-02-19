\name{osnd_fit}
\alias{osnd_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
osnd_fit
}
\description{
osnd_fit identifies OSND points on both original waveforms and model-generated waveforms, such that the HED
model's ability to recapitulate important morphological features and fiducial points can be assessed.
}
\usage{
osnd_fit(bf = beat_final, ppg, gs = model2.GetSegment, r = model2.Rebuild2, sf = splinefun, dp = diast_pk, oa = osnd_of_average, sr = samplingRate, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bf}{
finalised model parameter outputs
}
  \item{ppg}{
ppg time series
}
  \item{gs}{
model2.GetSegment function
}
  \item{r}{
model2.Rebuild2 function
}
  \item{sf}{
splinefun function
}
  \item{dp}{
diast_pk function
}
  \item{oa}{
osnd_of_average function
}
  \item{sr}{
sampling rate
}
  \item{plot}{
logical, plots each waveform with modelled version superimposed, as well as OSND points for each
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{osnd_diff }{list of values representing the error between model-generated OSND values and actual OSND values}
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
osnd_fits <- osnd_fit(beat_final, ppg, plot = F)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
