\name{GGplotFits}
\alias{GGplotFits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GGplotFits
}
\description{
GGplotFits produces GGplot plots of each modelled waveform in the batch, including component waves and baselines.
}
\usage{
GGplotFits(beats_in, ppg, beat2, gs = model2.GetSegment, rb = model2.Rebuild2, run, pr, p = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beats_in}{
number of beats in the batch
}
  \item{ppg}{
ppg time series
}
  \item{beat2}{
finalised dataframe of model parameter outputs for each waveform in the batch
}
  \item{gs}{
model2.GetSegment function
}
  \item{rb}{
model2.Rebuild2 function
}
  \item{run}{
participant number
}
  \item{pr}{
dose level
}
  \item{p}{
logical, setting to true stores the plot as a variable
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Automated plotting of each fitted waveform in the batch
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
GGplotFits(beats_in, ppg, beat2, gs = model2.GetSegment, rb = model2.Rebuild2, run, pr, p = F)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
