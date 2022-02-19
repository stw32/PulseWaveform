\name{sep_beats}
\alias{sep_beats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sep_beats
}
\description{
sep_beats segments the PPG time series into individual pulse waveforms, allowing for an average morphology to be elucidated as well as the later extraction of individual waveform features. Segmented beats undergo a cleaning process to exclude beats with aberrant morphology (thresholds have been determined empirically). There is also an optional subsetting routine which allows for identification of periods of increased heart rate (or lower IBI) in a time series and subsequent selection of the corresponding waveforms - for more details on this see supplementary material.
}
\usage{
sep_beats(odiff, bc, samp, wuv, wvlen, inx, o, ibi, scale = TRUE, q = FALSE, subset = FALSE, boundaries = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{odiff}{
vector of intervals between successive O (trough) points
}
  \item{bc}{
the baseline corrected time series
}
  \item{samp}{
sampling rate
}
  \item{wuv}{
dataframe of w, u and v values for each identified beat
}
  \item{wvlen}{
the median of O-diff +/- empirically determined value
}
  \item{inx}{
vector of inflection points (x/time coordinates)
}
  \item{o}{
vector of o points
}
  \item{ibi}{
vector of interbeat intervals (as determined by w points (peaks in 1st derivative))
}
  \item{scale}{
logical, determines if pulse waveforms are scaled / normalized for amplitude
}
  \item{q}{
logical, if TRUE number of beats excluded for each reason will be reported (useful for determining when thresholds may be to liberal/conservative)
}
  \item{subset}{
logical, determines if subsetting (as described above) is carried out
}
  \item{boundaries}{
ISO analysis specific argument (for carrying over subsetting constants between time series)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{average_wave }{vector of values for the averaged waveform}
\item{pulse }{dataframe of all individual pulse waveforms in segmented form}
\item{wuv }{dataframe of w, u and v values (updated for excluded beats) for all waveforms}
\item{rejects }{list of vectors, each comrprised of waves excluded for a particular reason}
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
tmp <- sep_beats(odiff = oDiff, bc = baseCor, samp = samplingRate, wuv = wuv, wvlen = waveLen,
                 ibi=ibi, o=o_orig, inx = inflexX, scale = T, q = F, subset = FALSE, boundaries)
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
