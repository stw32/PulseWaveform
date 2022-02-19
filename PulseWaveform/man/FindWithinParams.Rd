\name{FindWithinParams}
\alias{FindWithinParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FindWithinParams
}
\description{
FindWithinParams generates a 13x12 matrix for each waveform in an inputted batch of waveforms. Each row is a
unique parameter set for the waveform in question, while columns 1:12 correspond to each of the 12 parameters
(baseline1, baseline2, systolic time, systolic amplitude, systolic width, 2nd reflectance wave time, 2nd
reflectance wave amplitude, 2nd reflectance wave width, 1st reflectance wave time, 1st reflectance wave amplitude,
1st reflectance wave width, and decay rate).

Each populated row differs from the top row by a single parameter. The specific value that differs in each row
represents the optimised value for that parameter when all other parameters are held constant. Rows corresponding
to changes in parameters which are to be fixed across beats are unpopulated (NA). The matrices outputted by
FindWithinParams are the building blocks for the final matrix that will be inputted into the downhill simplex routine.
}
\usage{
FindWithinParams(beats_in, ppg, beat, gs = model2.GetSegment, fp = model2.FixParams3, ms = simplex.MakeSimplex3, m2 = model2.ChiSq3, beat_vector = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beats_in}{
number of beats in a batch, over which the width, reflectance timing, and decay rate parameters will be fixed
}
  \item{ppg}{
the preprocessed ppg time series
}
  \item{beat}{
the dataframe of model parameter outputs populated with initial starting parameters
}
  \item{gs}{
model2.GetSegment function
}
  \item{fp}{
model2.FixParams3 function
}
  \item{ms}{
simplex.MakeSimplex3 function
}
  \item{m2}{
model2.ChiSq3 function
}
  \item{beat_vector}{
a list consisting of 1. the number of beats in the batch 2. the starting points of each beat on the ppg time series 3. the end points of each beat on the ppg time series
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
\item{a }{list of n 13x12 matrices, where n = number of beats in the inputted batch}
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
within_params <- FindWithinParams(beats_in, ppg, beat = new_beat, gs = model2.GetSegment,
                                                fp = model2.FixParams3, ms = simplex.MakeSimplex3,
                                                m2 = model2.ChiSq3, beat_vector = beat_vector,
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
