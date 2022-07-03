\name{ArrangeOutputs}
\alias{ArrangeOutputs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ArrangeOutputs
}
\description{
ArrangeOutputs organises lists of outputs from various parts of the processing pipeline into a single output list.
}
\usage{
ArrangeOutputs(beat_final, beat_orig, features, pulse, fit_check, ps, pr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beat_final}{
dataframe of finalised model parameter outputs
}
  \item{beat_orig}{
dataframe of intially estimated model parameter outputs
}
  \item{features}{
dataframe of morphological features for each waveform
}
  \item{pulse}{
dataframe of all individual waveforms (discrete form)
}
  \item{fit_check}{
list of goodness of fit measures for each modelled waveform (ChiSq, Max error, NRMSE, aNRMSE)
}
  \item{ps}{
time series number (for ISO data)
}
  \item{pr}{
dose level (for ISO data)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{tmp }{list composed of:}
\item{beat_final }{as above}
\item{features }{as above}
\item{fit_check }{as above}
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
temp <- ArrangeOutputs(beat_final, beat_orig, features, pulse, fit_check, ps, pr)
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
