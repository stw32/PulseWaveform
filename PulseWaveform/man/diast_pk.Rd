\name{diast_pk}
\alias{diast_pk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
diast_pk
}
\description{
diast_pk finds the position of the diastolic peak (D) on the average wave so as to inform identification of OSND points on all individual waveforms. Since the average waveform being evaluated is also free from NA values, the resultant shift along the x-axis in terms of index points is also calculated.
}
\usage{
diast_pk(avw, sr, scale = F, dias_param = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{avw}{
average waveform
}
  \item{sr}{
sampling rate
}
  \item{scale}{
logical - to establish bounds of y-axis
}
  \item{dias_param}{
optional - to inform search for D if the waveform is a HED modelled one
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
\item{diastPk }{diastolic peak position on the average waveform}
\item{xShift }{difference in x-axis values of the average waveform due to removal of NA values}
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
tmp <- diast_pk(avw = avWave, sr = samplingRate, scale = T)
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
