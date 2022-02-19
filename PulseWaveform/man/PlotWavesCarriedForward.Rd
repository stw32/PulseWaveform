\name{PlotWavesCarriedForward}
\alias{PlotWavesCarriedForward}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PlotWavesCarriedForward
}
\description{
PlotWavesCarriedForward is a function specific to the ISO study data. The number of waves subsetted in the sep_beats function
is plotted for each participant (for each dose level) to allow for assessment of which samples are likely to be more or less
robust due to the number of beats included.
}
\usage{
PlotWavesCarriedForward(waves_carried_forward1, waves_carried_forward3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{waves_carried_forward1}{
list of number of beats subsetted for further analysis in the script, at iso dose level
}
  \item{waves_carried_forward3}{
list of number of beats subsetted for further analysis in the script, at saline dose level
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A histogram is automatically plotted indicating the distribution of frequency of subsetted beats across the dataset
%%  ~Describe the value returned
%%  If it is a LIST, use
%%\item{comp1 }{Description of 'comp1'}
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
PlotWavesCarriedForward(waves_carried_forward1, waves_carried_forward3)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
