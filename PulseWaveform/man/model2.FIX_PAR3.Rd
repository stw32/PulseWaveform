\name{model2.FIX_PAR3}
\alias{model2.FIX_PAR3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model2.FIX_PAR3
}
\description{
model2.FIX_PAR3 takes a set of 12 parameters as input and fixes any that exceed constraints defined within the function.
For every violation of model constraints a penalty is accrued, the sum of which is outputted along with the fixed parameters.
}
\usage{
model2.FIX_PAR3(time, within_beat_params, across_beat_params, debug=FALSE, renal_param, sys_t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
x values of the ppg segment being modelled
}
  \item{within_beat_params}{
parameters not fixed across beats
}
  \item{across_beat_params}{
parameters fixed across beats
}
  \item{debug}{
logical, if true function prints parameter values relative to constraint values, and the penalties accrued
}
  \item{renal_param}{
the timing of the 2nd reflectance wave, used as a reference to the initial estimation for this parameter
}
  \item{sys_t}{
the timing of the systolic wave, used as a reference to the initial estimation for this parameter
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{result }{vector of fixed parameters and total penalty accrued}
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
model2.FIX_PAR3 can handle parameter inputs that indicate only a single systolic wave, or two component waves only. Likewise
it can handle one rather than two baseline parameter inputs. This is to support modelling of waveforms with only as many
parameters as are necessary (hence also the use of reduced ChiSq as the chosen measure of goodness of fit). The optionality
of parameter inputs is summarized as follows:

params: {Baseline, {baseline 2}, t_sys, H_sys, W_sys, {dt_1, H_1, W_1, {dt_2, H_2, W_2}}}
across_beat_params: { w[1], t[2], w[2], t[3], w[3] }

Or:
par: {base1, {base2}, t[1], h[1], #, #, { h[2], #, #, { h[3], ..., ... }}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
temp <- model2.FIX_PAR3(time = data[, 1],
                          within_beat_params = params[c(1:4, 7, 10)],
                          across_beat_params, debug = F, renal_param = rp,
                          sys_t)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
