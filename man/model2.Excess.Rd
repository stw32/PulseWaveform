\name{model2.Excess}
\alias{model2.Excess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model2.Excess
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
model2.Excess(y, offset, baseline)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the y values of the waveform, as a vector
}
  \item{offset}{
the y value of the sample in the time series directly preceding the first value of the inputted waveform
}
  \item{baseline}{
anticpated baseline towards which the decay element decays
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{result }{the approximated excess element of the waveform in question}
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
residue <- e(data[,2], ppg[seg[1]-1,2], -0)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
