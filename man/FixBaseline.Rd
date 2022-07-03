\name{FixBaseline}
\alias{FixBaseline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FixBaseline
}
\description{
FixBaseline assesses the first and second fitted baselines for each waveform in a batch. If the baselines are close in value
and the reduced ChiSq value (goodness of fit) for a waveform's parameter set is lower when the baselines are equal (i.e
baseline2 equal to baseline1), then baseline2 will be made equal to baseline1.
}
\usage{
FixBaseline(new_beat, f = model2.ChiSq3, renal_param, dias_param, sys_time, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{new_beat}{
dataframe of model parameter outputs
}
  \item{f}{
model2.ChiSq3 function
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
\item{new_beat }{dataframe of model parameter outputs with baseline parameters further corrected}
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
new_beat <- FixBaseline(new_beat, f = model2.ChiSq4, renal_param, dias_param, sys_time, w)
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
