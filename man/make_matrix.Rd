\name{make_matrix}
\alias{make_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make_matrix
}
\description{
make_matrix combines individual waveform matrices for non-fixed parameters and the batch-wide matrix for fixed
parameters into a singular matrix for inputting into the downhill simplex routine. The dimensions of the matrix
are dependent on the number of beats in the inputted batch, such that matrix dimensions are:

   [(6*n + 6) + 1]   x    [(6*n + 6)]

Where n = number of beats in the inputted batch. Therefore a batch of 10 waveforms would generate a matrix of
67x66. The number of columns also represents the number of dimensions the simplex will navigate in optimising
parameters across the entire batch.
}
\usage{
make_matrix(sim, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim}{
the matrix representing variations in fixed parameters across the batch
}
  \item{a}{
the set of matrices representing variations in non-fixed parameters for each waveform in the batch
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{sim }{a combined matrix representative of both fixed and non-fixed parameters across the batch}
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
mat <- make_matrix(across_params, within_params)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
