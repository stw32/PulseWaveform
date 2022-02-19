\name{feature_extract}
\alias{feature_extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
feature_extract
}
\description{
feature_extract extracts descriptive morphological features for each waveform inputted.
}
\usage{
feature_extract(oa, p, pw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oa}{
a list of osnd values for each wave (x and y coordinates)
}
  \item{p}{
dataframe of all individual waveforms (discrete form)
}
  \item{pw}{
list of all individual waveforms (polynomial spline form)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{features }{a dataframe of all extracted features for each waveform}
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
features <- feature_extract(oa = osnd_all, p = pulse, pw = polyWave)
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
