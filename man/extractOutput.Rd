\name{extractOutput}
\alias{extractOutput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
extractOutput
}
\description{
extractOutput extracts the simplex-optimised values for fixed and non-fixed parameters across a given batch.
}
\usage{
extractOutput(beats_in, sim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beats_in}{
number of beats in the batch
}
  \item{sim}{
the final matrix outputted by the downhill simplex routine
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{temp }{list consisting of:}
\item{across }{optimised fixed parameters across the entire batch (6 in total)}
\item{within }{list of optimised non-fixed parameters for each waveform in the batch}
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
output <- extractOutput(beats_in, sim)
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
