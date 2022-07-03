\name{FactorAdjust}
\alias{FactorAdjust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FactorAdjust
}
\description{
FactorAdjust identifies the most appropriate factor value to be used in undetrending an ISO study ppg time
series (for inputting into the Undetrend function). It focuses on a single beat segment of the ppg time
series and adjusts the factor value to correct morphological distortions common to detrending in the ISO
data. These are 1. abnormal tail decay (e.g positive gradient to tail) 2. values significantly below the
approximated baseline (such that notch features are negative and therefore meaningless). The function
operates within the FindUndetrendingParams function.
}
\usage{
FactorAdjust(data, factorCutoff, ppg, u, beat, a., test, gs=gs, beatTime, nextTime, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  an individual beat segment of raw ppg time series data
}
  \item{factorCutoff}{
  the maximum acceptable tail gradient (default 0)
}
  \item{ppg}{
  the ppg time series
}
  \item{beat}{
  dataframe of all detected peaks in the ppg time series
}
  \item{u}{
  Undetrend function
}
  \item{a.}{
  index number indicating which beat (relative to the calculated minimum IBI of the time series) to sample
}
  \item{test}{
  index number indicating the beat identified as the minimum IBI of the time series
}
  \item{gs}{
  model2.GetSegment function
}
  \item{beatTime}{
  index value for the start of the beat to be assessed
}
  \item{nextTime}{
  index value for the end of the beat to be assesssed
}
  \item{plot}{
  logical, if TRUE, plots the beat segment for each iteration of the while loop, demonstrating the change
  in morphology as it is corrected
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{factor_value }{the appropriate factor value for restoring proper morphology to the ppg time series, to be inputted into the Undetrend function to complete ISO study data preprocessing}
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
 A minimum IBI section of the time series is relevant to this function given its purpose in preprocessing ISO
 study data. The morpological change associated with isoprenaline results in a lowered notch (and IBI), and thus
 any alterations to morphology to ensure notches above 0 should be targeted at the waves with the lowest notches.

 Within the ISO data four types of waveform endings are possible:

  1. prolonged negative slope (morphologically normal)
  2. prolonged positive slope (due to heavy detrending)
  3. ~5 positive slope coming out of the notch, followed by ~5 negative slope as a rather short tail i.e an n shape
  (these tend to be ISO waves (they are short due to high HR))
  4. ~5 positive slope at the end of the tail, preceded by a prolonged negative slope i.e a v shape (due to tail noise)

 1. and 2. will be processed the same regardless of which of the last 10 datapoints the gradient is taken from
 3. In the initial 'last 5' datapoint assessment, these will be found to be negative and thus will continue being
 iterated on based on the last 5 datapoints.
 4. In the initial 'last 5' datapoint assessement, these will be found to be positive, and so will be iterated on
 based on NOT the last 5 datapoints.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
FactorAdjust(data, factorCutoff, ppg, u, beat, a., test, gs=gs, beatTime, nextTime, plot = T)
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
