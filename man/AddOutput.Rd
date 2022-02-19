\name{AddOutput}
\alias{AddOutput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AddOutput
}
\description{
AddOutput adds output columns to the beat dataframe (see general purpose script) to be populated with model
parameter outputs.
}
\usage{
AddOutput(beat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beat}{
dataframe where number of rows = number of waveforms in sample
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{beat }{beat dataframe with output columns added}
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
beat <- AddOutput(beat)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
