\name{model2.Peak}
\alias{model2.Peak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model2.Peak
}
\description{
model2.Peak generates a sine wave from the three parameters used to define each component wave (timing, width,
and amplitude).
}
\usage{
model2.Peak(time, peakParams)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
x-axis values of the sine wave
}
  \item{peakParams}{
vector of three parameters defining a given peak
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{result }{sine wave constructed from model peak parameters}
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
result <- result + model2.Peak(xy[,1],params[3:5])
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
