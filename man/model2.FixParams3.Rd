\name{model2.FixParams3}
\alias{model2.FixParams3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model2.FixParams3
}
\description{
model2.FixParams3 organises inputs and outputs of the model2.FIX_PAR3 function.
}
\usage{
model2.FixParams3(data,params, across_beat_params = NULL, debug=FALSE, rp = renal_param, sys_t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
ppg time series segment relevant to waveform being modelled
}
  \item{params}{
vector of model parameters
}
  \item{across_beat_params}{
vector of fixed parameters (set to null if all parameters to be sourced from params input)
}
  \item{debug}{
logical, now redundant
}
  \item{rp}{
1st reflectance wave timing paramter intiial estimation
}
  \item{sys_t}{
systolic wave timing paramter intiial estimation
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{temp }{vector of fixed parameters}
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
fixed <- FixOutput(beats_in, beat = new_beat, ppg, gs = model2.GetSegment,
                                 fp = model2.FixParams3, across = output[[1]], within = output[[2]],
                                 sys_time = sys_time)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
