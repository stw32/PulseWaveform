\name{model2.ChiSq3}
\alias{model2.ChiSq3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model2.ChiSq3
}
\description{
model2.ChiSq3 calculates a reduced ChiSq (goodness of fit) value for a model generated waveform, when provided
with the original waveform being modelled and a set of model parameters. The function can either be applied to
calculate the goodness of fit of an individual waveform, or batches of waveforms, and is used to drive the
downhill simplex optimisation of model parameters.
}
\usage{
model2.ChiSq3(data, params, debug=FALSE, beats, optional = NULL, beat = NULL, a = NULL, plot = FALSE, renal_param, dias_param, sys_time, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
section of ppg time series
}
  \item{params}{
model parameters
}
  \item{debug}{
logical, currently redundant
}
  \item{beats}{
list of number of inputted beats, there beginnings in the time series, and their endings
}
  \item{optional}{
logical, now redundant
}
  \item{beat}{
dataframe of model parameters
}
  \item{a}{
the combined matrix used in the downhill simplex routine (alternative source of parameters)
}
  \item{plot}{
logical, if set to true plots model generated waveform against original waveform
}
  \item{renal_param}{
the starting parameter for 1st reflectance peak timing (inputted to prevent drastic deviations from this value)
}
  \item{dias_param}{
the starting parameter for 2nd reflectance peak timing (inputted to prevent drastic deviations from this value)
}
  \item{sys_time}{
the starting parameter for systolic peak timing (inputted to prevent drastic deviations from this value)
}
  \item{w}{
the timing of the 1st derivative peaks on the ppg time series
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{ts_fit }{total reduced ChiSq value summed across all beats inputted}
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
makesimplex3 inputs a single set of parameters (tParam) to test, whereas run.simplex2 inputs a matrix
This determines where within beat parameters should be extracted from, hence beat and a are NULL unless otherwise specified.

During the downhill simplex parameter optimization process, movements of parameter values outside of defined contraints
will lead to the application of penalties, which model2.ChiSq3 adds to the ChiSq value outputted, driving the
optimisation process to occur within constraints. The penalties are proportionate to the degree to which the constraints
are exceeded, to prevent mere ignorance of the constraint once violated.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
chiSq[i] <- f(data, params = NULL, optional=NULL,
                  a = result[i, ], beats = beat_vector,
                  renal_param = renal_param,
                  dias_param = dias_param,
                  sys_time = sys_time, w = w)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
