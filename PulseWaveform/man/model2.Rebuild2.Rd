\name{model2.Rebuild2}
\alias{model2.Rebuild2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model2.Rebuild2
}
\description{
model2.Rebuild2 takes HED model parameters and reconstructs a segment of the ppg time series corresponding to one waveform.
First the excess is calculated from the 9 peak parameters, then the decay element is added to generate a final waveform.
}
\usage{
model2.Rebuild2(xy,offset,params,invert=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{
x and y values for a particular segment of the ppg time series
}
  \item{offset}{
the y value of the data point immediately preceding the inputted segment of the ppg time series
}
  \item{params}{
vector of 12 parameters for generating a modelled wave
}
  \item{invert}{
logical, if set to true adds decay element as well as excess element
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{result }{dataframe of x and y values defining the model generated waveform}
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
PlotFits(beats_in, ppg, beat2, gs = model2.GetSegment, rb = model2.Rebuild2)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
