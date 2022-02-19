\name{GetDirec}
\alias{GetDirec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GetDirec
}
\description{
GetDirec finds the file directory for an individual participant
}
\usage{
GetDirec(run, Participants, dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{run}{
participant number
}
  \item{Participants}{
list of participants
}
  \item{dir}{
file directory for the ISO dataset
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{direc }{file directory for an individual participant}
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
direc <- GetDirec(run, Participants, dir)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
