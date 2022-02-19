\name{find_sd}
\alias{find_sd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
find_sd
}
\description{
find_sd is a clone of find_average; only the output is changed (from mean to SD)
}
\usage{
find_sd(p, ao)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
dataframe consisting of all segmented waveforms
}
  \item{ao}{
a vector indicating the end point of each waveform
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{sd_wave }{a vector consisting of standard deviation values (can also be plotted to indicate variance)}
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
sd_wave <- find_sd(p = pulse, ao = afterO)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
