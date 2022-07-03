\name{PlotRejects}
\alias{PlotRejects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PlotRejects
}
\description{
PlotRejects takes a list of the output list 'reject' from the sep_beats function and consutructs an R plot indicating the number
of rejects within a given samples, colour coded to indicate specific reasons for rejection. The function can be used as a tool
to assess how conservative the sep_beats function is behaving for a given sample such that appropriate adjustments to
empirical thresholds can be altered if necessary. Initially written to work within the ISO study main script.
}
\usage{
PlotRejects(rejected_waves_list1, rejected_waves_list3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rejected_waves_list1}{
a list of number of beats rejected for each reason, for a given sample, nested within a list of samples (participants in ISO study)
}
  \item{rejected_waves_list3}{
equivalent data structure to rejected_waves_list1 but for a different category of sample (e.g dose level in ISO study) (optional)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Automatically plots when called.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
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
The code for this function is inefficient and, while functioning, will need cleaning up for future releases.

Reasons for beat rejection include 1. excessively long waveforms (dotted), 2. excessively short waveforms (red), 3. bi-peaked
waveforms (blue), 4. waveforms that begin to enter a second systolic upstroke (green), 5. waveforms that drop significantly
below baseline (orange), 6. waveforms with significant variability in relation to the average (brown), 7. waveforms that are
significantly different in morphology than the average (purple), and total number of rejected beats (black))
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
PlotRejects(rejected_waves_list1, rejected_waves_list3)

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
