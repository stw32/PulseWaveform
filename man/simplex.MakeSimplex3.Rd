\name{simplex.MakeSimplex3}
\alias{simplex.MakeSimplex3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simplex.MakeSimplex3
}
\description{
simplex.MakeSimplex3 iterates on the initially estimated values for non-fixed parameters, refining them to give the simplex an improved starting position. It is used within the FindWithinParams function, such that a matrix of optimised non-fixed parameters can be generated for each beat in an inputted batch.
}
\usage{
simplex.MakeSimplex3(ppg, param,f,inScale, directions=NULL, inTol=-1, optional=NULL, debug=FALSE, beat_vector = beat_vector, renal_param, dias_param = dias_param, sys_time, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppg}{
ppg time series segment
}
  \item{param}{
model parameter inputs
}
  \item{f}{
model2.ChiSq3 function
}
  \item{inScale}{
the starting value by which to iterate on parameter values (0.1 recommended)
}
  \item{directions}{
the initial direction in which to iterate parameter values (default null)
}
  \item{inTol}{
the minimum value a parameter can change if it cannot be iterated upon to improve goodness of fit
}
  \item{optional}{
for direct input into model2.ChiSq3 function, default NULL
}
  \item{debug}{
logical, if true iteration and refinement statuses will be printed to console
}
  \item{beat_vector}{
an index of beats to be modeled and their x-coordinates to extract from the PPG time series (for direct input to f)
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
\item{result }{a matrix of parameter sets, including the initially inputted parameters (row 1) and a parameter set for each refined parameter}
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
simplex.MakeSimplex3 is for non-fixed (within beat) parameters only. Thus the outputted matrix will include absent (NA) rows where non-fixed parameters are represented (these are refined separately in the simplex.MakeSimplex2 function).
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
