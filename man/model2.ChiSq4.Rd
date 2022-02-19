\name{model2.ChiSq4}
\alias{model2.ChiSq4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model2.ChiSq4
}
\description{
model2.ChiSq4 is a clone of model2.ChiSq3. It differs by providing additional measures of goodness of fit besides reduced
ChiSq for the inputted batch. These are listed below.
}
\usage{
model2.ChiSq4(data, params, debug=FALSE, beats, beat, a = NULL, plot = FALSE, renal_param, dias_param, sys_time, w)
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
\item{fit }{list consisting of:}
\item{ts_fit }{total reduced ChiSq value summed across all beats inputted}
\item{beat_fit }{list of reduced ChiSq values for each individual beat}
\item{max_error }{residual of greatest value for each individual beat}
\item{NRMSE }{Normalised root mean square error (see supplementary material)}
\item{aNRMSE }{alternative normalised root mean square error (see supplementary material)}
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
fit_check[[k]] <- model2.ChiSq4(data = ppg, params = NULL, beats = beat_vector,
                                            beat = new_beat, a = sim[1, ], plot = FALSE,
                                            renal_param = renal_param, dias_param = dias_param,
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
