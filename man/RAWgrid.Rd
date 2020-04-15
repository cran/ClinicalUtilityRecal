\name{RAWgrid}
\alias{RAWgrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RAW grid for Cross Validation
}
\description{
Tuning parameters for calibration weights may not be intutitive. Instead, relative average weight (RAW) can be used to find sensible tuning parameters. The RAW is the average weight of observations within the RAW region, \eqn{[R_{l,RAW},R_{u,RAW}]}, divided by the average weight of observations outside RAW region.

This function elicits a grid of tuning parameter \eqn{lambda} or \eqn{delta} given a sequence of relative average weights (RAW) values, to be used for cross-validation. One tuning parameter \eqn{delta} or \eqn{lambda} must be specified.

Two clinically relevant lower and upper bound values must be specified. The first, \eqn{[R_{l},R_{u}}, define the clinically relevant region used for the weighting function. The second, \eqn{[R_{l,RAW},R_{u,RAW}]}, define the complemenatry regions used for defining relative average weights.
}
\usage{
RAWgrid(r, rl, ru, p, y, rawSeq=seq(0.1,0.9,0.1), cvParm,
delta = NULL, lambda = NULL, rl.raw, ru.raw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{r}{
  Clinically relevant risk threshold}
  \item{rl}{
  Lower bound of clinically relevant region}
    \item{ru}{
  Upper bound of clinically relevant region}
    \item{p}{
  Vector of risk score values}
    \item{y}{
  Vector of binary outcomes, with 1 indicating event (cases) and 0 indicating no event (controls)}
  \item{rawSeq}{
  Sequence of relative average weights}
  \item{cvParm}{
  Tuning parameter that will be selected via cross-validation. Can either be \eqn{lambda} or \eqn{delta}
}
  \item{delta}{
  Calibration weight for observations outside clinically relevant region \eqn{[R_l,R_u]}. Needs to be specified if cvParm="lambda", otherwise should be left blank, and function will generate a sequence of delta for cross-validaton.
}
  \item{lambda}{
Tuning parameter for controlling exponential decay of calibration weights. Needs to be specified if cvParm="delta", otherwise should be left blank, and function will generate a sequence of lambda for cross-validaton.}
  \item{rl.raw}{
Lower bound for defining the relative average weight region of interest
}
  \item{ru.raw}{
Upper bound for defining the relative average weight region of interest}
}

\value{
\item{rrWt.seq }{Matrix containing the input RAW sequence, and corresponding weight tuning parameter lambda and weight value delta (with one fixed depending on cvParm selection). If 'NA' is returned there may be too few events within the weight clinically relevant RAW region, meaning that RAW value is not possible}
}
\references{Mishra, A. (2019). Methods for Risk Markers that Incorporate Clinical Utility (Doctoral dissertation). (Available Upon Request)}



\author{
Anu Mishra
}
\section{Warning }{
{Event rate outside RAW interval, widen RAW interval}{If no or too cases are inside RAW interval [R_{l},R_{u}], cross-validation procedure may not be stable, so this warning indicates RAW interval should widen to include more cases.}
}


\examples{
### Load data ##
data(fakeData)

### Get Grid of Tuning Parameters ###
grid <- RAWgrid(r = 0.3,rl = -Inf,ru = Inf,p = fakeData$p,y = fakeData$y,
                cvParm = "lambda",delta=1,rl.raw = 0.25,ru.raw = 0.35)
}


\seealso{
\code{\link{calWt}},
\code{\link{cvWtTuning}},
\code{\link{cvRepWtTuning}}
}
