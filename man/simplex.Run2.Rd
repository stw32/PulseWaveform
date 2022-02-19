\name{simplex.Run2}
\alias{simplex.Run2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simplex.Run2
}
\description{
Simplex.Run2 intiates the Nelder Mead downhill simplex routine to optimise model parameters and goodness of fit.
}
\usage{
simplex.Run2(data = ppg,simplexParam = mat, f = model2.ChiSq3, optional=NULL, beat_vector = beat_vector, ms = simplex_iterations, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w, run = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
ppg time series segment
}
  \item{simplexParam}{
the matrix outputted by make.matrix, representative of all model parameters in an inputted batch of beats
}
  \item{f}{
model2.ChiSq3 function
}
  \item{optional}{
currently redundant
}
  \item{beat_vector}{
an index of beats to be modeled and their x-coordinates to extract from the PPG time series
}
  \item{ms}{
number of simplex iterations to conduct
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
  \item{run}{
indicates which iteration of the simplex routine is being run
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix with equivalent structure to the one inputted, with optimised parameter values.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
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
sim <- simplex.Run2(data = ppg, simplexParam = mat, f = model2.ChiSq3, optional=NULL,
                                  beat_vector = beat_vector, renal_param = renal_param,
                                  dias_param = dias_param, sys_time = sys_time, w = w, run = c("run", i))
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
