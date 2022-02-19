\name{UpdateBeat}
\alias{UpdateBeat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
UpdateBeat
}
\description{
UpdateBeat generates a new matrix of model parameters by replacing values in the existing dataframe with
parameters that are both optimised and constrained.
}
\usage{
UpdateBeat(beats_in, beat, fixed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beats_in}{
number of beats in the batch
}
  \item{beat}{
existing dataframe of model parameter estimates (populated with initial estimations)
}
  \item{fixed}{
list of optimised and constrained parameters for each waveform, to replace initial estimates
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{new_beat }{Dupdated dataframe of model parameter outputs}
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
new_beat <- UpdateBeat(beats_in, beat, fixed)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
