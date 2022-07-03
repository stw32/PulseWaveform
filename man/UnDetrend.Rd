\name{UnDetrend}
\alias{UnDetrend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
UnDetrend
}
\description{
UnDetrend is a function used as part of ISO data preprocessing. It reverses a detrending effect apparent when
viewing the data and presumably originating from the hardware source, according to the inputted factor and
offset values.
}
\usage{
UnDetrend(ppg,factor,offset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppg}{
ppg time series to be pre-processed  / undetrended
}
  \item{factor}{
determines the degree of y-axis decay between data points in the time series
}
  \item{offset}{
determines the overall trend of the time series
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{result }{an undetrended ppg time series}
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
ppg[,2] = UnDetrend(ppg,factor=factor_value,offset=offset_value)
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
