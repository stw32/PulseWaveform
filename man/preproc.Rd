\name{preproc}
\alias{preproc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
preproc
}
\description{
preproc undertakes a pre-processing routine specific to raw Bioradio data. This entails downsampling the data and then 'undetrending' it. Downsampling removes repeated values (the BioRadio device provides 250 samples per second, but the PPG is only sampled 75 times per second), whilst undetrending reverses a detrending function inherent to Bioradio hardware. To expand on the latter, analysis of device output indicates that the PPG signal is detrended by application of the following formula:

              OUT[i] = 80 + (OUT[i-1]-80) * 0.96875 + (IN[i] - [IN[i-1])

where the constant 0.96875 is an approximation fitted to the data. Individual pulse events are more comprehensible once undetrended.
}
\usage{
preproc(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
raw Bioradio data as outputted by device
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A downsampled and undetrended time series
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
data <- read.csv(file.choose(), header = T)
undetrended_data <- data.frame(preproc(dat=data))
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
