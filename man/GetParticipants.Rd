\name{GetParticipants}
\alias{GetParticipants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GetParticipants
}
\description{
GetParticipants creates a list of all useable files in the ISO dataset
}
\usage{
GetParticipants(direc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{direc}{
file directory for the ISO dataset
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{string_list }{list of ISO participant IDs}
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
Participants <- GetParticipants(dir)
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
