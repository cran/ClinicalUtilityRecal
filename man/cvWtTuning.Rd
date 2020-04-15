\name{cvWtTuning}
\alias{cvWtTuning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation for Selecting Weight Tuning Parameter
}
\description{
Calibration weights require specification of tuning parameter \eqn{delta} or \eqn{lambda}. This function uses K-fold cross-validation to select tuning parameter used for calibration weights, with standardized net benfeit (sNB) as objective function. Either one of \eqn{delta} or \eqn{lambda} must be specificed. The sequence of tuning parameters can be obtained from the \code{\link{RAWgrid}} function.
}
\usage{
cvWtTuning(p,y,r,rl,ru,kFold=5,cvParm,tuneSeq,cv.seed=1111)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{
  Vector of binary outcomes, with 1 indicating event (cases) and 0 indicating no event (controls)}
    \item{p}{
  Vector of risk score values}
    \item{r}{
  Clinically relevant risk threshold}
  \item{rl}{
  Lower bound of clinically relevant region}
    \item{ru}{
  Upper bound of clinically relevant region}
    \item{kFold}{
  Number of folds for cross-validation}
  \item{cvParm}{
  Parameter to be selected via cross-validation. Can be either \eqn{delta} the weight assigned to observations outside the clinically relevant region [R_{l},R_{u}], or the \eqn{lambda} tuning parameter controlling exponential decay within the clinically relevant region [R_{l},R_{u}]}
  \item{tuneSeq}{
  Sequence of values of tuning parameters to perform cross-validation over}
  \item{cv.seed}{
  Intial seed set for random splitting of data into K folds}

}
\value{
\item{cv.res}{Matrix containing sequence of tuning parameters and corresponding cross-validation sNB}
\item{cv.param}{Value of tuning parameter selected via cross validation}
\item{cv.full}{Matrix of cross-validation results for all folds}
}

\references{Mishra, A. (2019). Methods for Risk Markers that Incorporate Clinical Utility (Doctoral dissertation). (Available Upon Request)}



\author{
Anu Mishra
}
\note{
Note this function does not split data into training and validaion set, but performs the K-fold cross-validation procedure on all data included. We advise that a separate, validation subset should be split from the data used in this function.
}

\seealso{
\code{\link{calWt}},
\code{\link{RAWgrid}},
\code{\link{nb}},
\code{\link{cvRepWtTuning}}

}
