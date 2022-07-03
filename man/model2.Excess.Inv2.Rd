\name{model2.Excess.Inv2}
\alias{model2.Excess.Inv2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model2.Excess.Inv2
}
\description{
model2.Excess.Inv2 adds the model's decay element to its excess element. It is used within model2.Rebuild2 as part of the
construction of a modelled waveforms from inputted model parameters.
}
\usage{
model2.Excess.Inv2(time,excess,offset,baselineStart,baselineEnd,timeBase,config.rate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
x coordinates of a particular segment of the PPG time series
}
  \item{excess}{
segment corresponding to excess element of reconstructed waveform
}
  \item{offset}{
the y value of the data point immediately preceding the inputted segment of the ppg time series
}
  \item{baselineStart}{
the y-value towards which the first baseline decays
}
  \item{baselineEnd}{
the y-value towards which the second baseline decays
}
  \item{timeBase}{
the point in time (x value) when baseline 1 switches to baseline 2, defined as the D (2nd reflectance) peak timing
}
  \item{config.rate}{
decay rate parameter
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{result }{a modelled waveform with combined excess and decay elements}
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
result <- model2.Excess.Inv2(xy[,1],result,offset,params[1],params[2],
                                 params[3]+1*params[6],
                                 config.rate = params[12])
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
