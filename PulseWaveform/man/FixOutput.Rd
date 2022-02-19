\name{FixOutput}
\alias{FixOutput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FixOutput
}
\description{
FixOutput applies constraints to the optimised model parameter outputs to ensure they remain within feasible
limits, by applying the model2.FixParams3 function.
}
\usage{
FixOutput(beats_in, beat, ppg, gs = model2.GetSegment, fp = model2.FixParams3, across = output[1], within = output[2], sys_time = sys_time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beats_in}{
number of beats in an inputted batch
}
  \item{beat}{
dataframe of model output parameters (not yet updated with optimised parameters)
}  \item{ppg}{
ppg time series (unsegmented)
}  \item{gs}{
model2.GetSegment function
}  \item{fp}{
model2.FixParams3 function
}  \item{across}{
optimised fixed parameters across the batch
}  \item{within}{
optimised non-fixed parameters for each waveform in the batch
}  \item{sys_time}{
the starting parameter for systolic peak timing (inputted to prevent drastic deviations from this value)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{fixed }{list of vectors, each vector corresponds to the 12 optimised and fixed parameters for one waveform in the batch}
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
