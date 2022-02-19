\name{preclean_wuv}
\alias{preclean_wuv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
preclean_wuv
}
\description{
Preclean_wuv takes u, v and w values and assesses the x-axis (time) difference between them for each beat. For normal beats, time from u to w is around half (50 percent) of the time from u to v. Beats with abnormal / artefactual systolic upstrokes tend to have outlying values for this measure. Thus PreClean can identify and remove them.
}
\usage{
preclean_wuv(w, uv, o, samp, sp, q = FALSE)
}
\arguments{
  \item{w}{
vector of all x-axis coordinates corresponding to w points
}
  \item{uv}{
a dataframe containing x and y coordinates for all u and v values
}
  \item{o}{
vector of o points
}
  \item{samp}{
sample rate
}
  \item{sp}{
The original PPG time series, discrete form
}
  \item{q}{
Logical, will pause function and give the option of plotting rejected beats
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{w }{as above, with aberrant beats removed}
\item{uv }{  "       "     }
\item{o }{  "       "     }
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
tmp <- preclean_wuv(w=w, uv=uv, o=o_orig, samp = samplingRate, sp = spline1, q = F)
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
