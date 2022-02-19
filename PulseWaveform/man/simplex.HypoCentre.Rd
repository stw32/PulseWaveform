\name{simplex.HypoCentre}
\alias{simplex.HypoCentre}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simplex.HypoCentre
}
\description{
simplex.HypoCentre identifies the average vertex of the simplex when the apex (highest point / parameter set yielding worst
goodness of fit) is excluded. This is referred to as the 'node' or 'hypocentre' and represents the point about which the
worst vertex is reflected in order to generate a new and lower vertex.
}
\usage{
simplex.HypoCentre(mat_Param, index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat_Param}{
a simplex matrix
}
  \item{index}{
the row of the inputted matrix corresponding to the worst (poorest fit) parameter set
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{result }{a vertex representing the mean of simplex vertices when the apex is excluded (equivalent to a triangle's base)}
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
node <- simplex.HypoCentre(result,high)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
