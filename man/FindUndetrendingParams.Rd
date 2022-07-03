\name{FindUndetrendingParams}
\alias{FindUndetrendingParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FindUndetrendingParams
}
\description{
FindUndetrendingParams is specific to the ISO study analysis and preprocesses the data therein. However, it
may also be useful for preprocessing any PPG time series data in need of morphological correction due to
hardware detrending algorithms. It functions by optimising factor and offset undetrending parameters across
a number of sample waveforms from the four PPG time series per participant in the ISO study. From this
reasonable values to apply to all four time series are inferred.
}
\usage{
FindUndetrendingParams(direc, gs = model2.GetSegment, u = UnDetrend, oa = OffsetAdjust, fa = FactorAdjust, factorCutoff = 0, sr = samplingRate, pk_thrshd, pairs, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{direc}{
the directory for an individual participant in the ISO study
}
  \item{gs}{
model2.GetSegment function
}
  \item{u}{
UnDetrend function
}
  \item{oa}{
OffsetAdjust function
}
  \item{fa}{
FactorAdjust function
}
  \item{factorCutoff}{
the maximum acceptable tail gradient (default 0)
}
  \item{sr}{
sample rate
}
  \item{pk_thrshd}{
the objective threshold for initial identification of peaks in the 1st derivative (see ISO main script)
}
  \item{pairs}{
file names of the two pairs of time series used for each participant in the ISO study analysis
}
  \item{plot}{
logical, if TRUE plots data during both factor and offset adjustment
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{values }{a vector containing the offset and factor values required to undetrend the raw PPG data}
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
temp <- FindUndetrendingParams(direc = direc[1], gs = model2.GetSegment, oa = OffsetAdjust,
                                     fa = FactorAdjust, u = UnDetrend, factorCutoff = 0, sr = samplingRate,
                                     pairs = pairs, pk_thrshd = pk_thrshd, plot = F)
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
