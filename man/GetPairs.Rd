\name{GetPairs}
\alias{GetPairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GetPairs
}
\description{
GetPairs is an ISO analysis specific function. It is used to identify and order the four data files per
participant in the ISO dataset, arranding them into two pairs of 'isoprenaline' time series and 'saline'
time series.
}
\usage{
GetPairs(direc, run_order, participant_number, subjectID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{direc}{
the directory for an individual participant in the ISO dataset
}
  \item{run_order}{
the order of infusions given during the ISO study protocol
}
  \item{participant_number}{
the individual participant number for a given participant
}
  \item{subjectID}{
the participant ID as coded in the ISO study data files
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{pairs }{list of two vectors, each vector containing one ISO time series and one Saline time series}
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
pairs <- GetPairs(direc = direc[1], run_order, participant_number = run, subjectID = direc[2])
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
