\name{FindStartParams}
\alias{FindStartParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FindStartParams
}
\description{
FindStartParams generates starting parameters for the HED model, using estimates derived from the data itself.
Since the HED model assumes both decay and excess elements, an approximation of the decay is removed from
the original waveform data to infer the excess. Once the excess is defined, peaks in the excess are
searched for within empirically defined ranges. Values for width, timing and amplitude for each peak can then
be inferred.
}
\usage{
FindStartParams(batch_number, beats_in, beat, ppg, gs = model2.GetSegment, e = model2.Excess, sep = model2.SubtractExcessPeak, o_points = inflexX[o_orig], wuv = wuv, inflexX = inflexX, all_beats = FALSE, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{batch_number}{
the number of batches of waves to process
}
  \item{beats_in}{
the number of beats in a batch, over which the width, reflectance timing, and decay rate parameters will be fixed
}
  \item{beat}{
empty dataframe to be filled with starting parameters
}
  \item{ppg}{
ppg time series (unsegmented)
}
  \item{gs}{
model2.GetSegment function
}
  \item{e}{
model2.Excess function
}
  \item{sep}{
model2.SubtractExcessPeak function
}
  \item{o_points}{
vector of locations of the o (origin) points for each waveform in the ppg time series
}
  \item{wuv}{
dataframe of locations of the w, u and v points for each waveform in the ppg time series
}
  \item{inflexX}{
vector of x coordinates for each inflection point on the ppg time series
}
  \item{all_beats}{
logical, if true ensures all available beats in the time series are modeled
}
  \item{plot}{
logical, if true plots each beat for which parameters are being estimated
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{temp }{list consisting of the following:}
\item{beat }{updated dataframe populated with estimated starting parameters, to be used as input for downhill simplex}
\item{ppg }{additional columns providing parallel time series of inferred baseline, excess and residuals of ppg data}
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
temp <- FindStartParams(batch_number, beats_in, beat, ppg, gs = model2.GetSegment,
                                  e = model2.Excess, sep = model2.SubtractExcessPeak,
                                  o_points = inflexX[o_orig], wuv = wuv, inflexX = inflexX, all_beats)
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
