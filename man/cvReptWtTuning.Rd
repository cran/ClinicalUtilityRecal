\name{cvRepWtTuning}
\alias{cvRepWtTuning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Repeated Cross Validation for Weight Tuning Parameter Selection
}
\description{
Calibration weights require specification of tuning parameter \eqn{delta} or \eqn{lambda}. Since a single round of cross-validation can be noisy, cross-validation can be repeated multiple times with independent random partitions and the results be averaged. This function implements a repeated K-fold cross-validation where tuning parameter \eqn{labmda} or \eqn{delta} is selected by maximizing standardized net benefit (sNB) (i.e. repeated \code{\link{cvWtTuning}} procedure).

A a "one-standard error" rule can be used for selecting tuning parameters. Under the “one-standard error" rule the calibration weight tuning parameter (\eqn{lambda} or \eqn{delta}) is selected such that corresponding cross-validated sNB is within one-standard deviation of the maximum cross-validated sNB. This provides protection against overfitting the data and selecting a tuning parameter that is too extreme. If the "one-standard error" rule is not implemented, then the tuning parameter with the larged average cross-validted sNB (across folds and repetition) will be selected.
}
\usage{
cvRepWtTuning(y,p,r,rl,ru,kFold=5,cvRep=25,cvParm,tuneSeq,stdErrRule=TRUE,int.seed=11111)
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
  \item{cvRep}{
  Number of cross-validation repititions}
  \item{cvParm}{
  Parameter to be selected via cross-validation. Can be either \eqn{delta} the weight assigned to observations outside the clinically relevant region [R_{l},R_{u}], or the \eqn{lambda} tuning parameter controlling exponential decay within the clinically relevant region [R_{l},R_{u}]}
  \item{tuneSeq}{
  Sequence of values of tuning parameters to perform cross-validation over}
  \item{stdErrRule}{
  Use "one-standard" error rule selecting tuning parameter}
  \item{int.seed}{
  Intial seed set for random splitting of data into K folds}
}
\details{
To estimate the standard deviation of the cross-validated sNV, the dependence between the different partitions of cross-validation needs to be accounted for. Gelman (1992) give a variance estimator of convergence diagnostic statistic used when Markov Chain Monte Carlo with multiple chains are performed. The variance estimator accounts for both the variability of the statistic “within" a single chain, and the variance of the statistic across, or “between", chains. Analogously, we can use this framework to estimate the “within" repetition variance (i.e. variation in sNB from a single round of K-fold cross-validation) and the “between" repetition variance. We denote the ‘within" repetition variance as W and the “between" repetition variance as B . We augment this formula slightly from that given in Gelman (1992) to account for the fact that as the number of cross-validation repetitions increases, the between-repetition variability should decrease. See Mishra et al (2020) for full expressions of B and W.}
\value{
\item{cv.sNB}{Standardized net benefit (sNB) of tuning parameter selected via cross-validatoin }
\item{cv.RAW}{Corresponding RAW value given cross-valiated selected tuning parameter}
\item{cv.lambda}{\eqn{lambda} value selected via cross-validation if \eqn{cvParm=lambda}, otherwise user specified \eqn{lambda} value}
\item{cv.delta}{\eqn{delta} value selected via cross-validation if \eqn{cvParm=delta}, otherwise user specified \eqn{lambda} value}
\item{avgCV.res}{Averaged (across-replications) cross-validated sNB for sequence of tuning parameters}
\item{W}{Estimate of "with-in" repetition variance. Will only return if stdErrRule==TRUE}
\item{B}{Estimate of "between" repetition variance. Will only return if stdErrRule==TRUE}

\item{fullList}{List of cross-valiation results for all fold and repititions}

}
\references{Mishra, A. (2019). Methods for Risk Markers that Incorporate Clinical Utility (Doctoral dissertation). (Available Upon Request)

Friedman, J., Hastie, T., & Tibshirani, R. (2001). The elements of statistical learning (Vol. 1, No. 10). New York: Springer series in statistics.

Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. Statistical science, 7(4), 457-472.
}
\author{
Anu Mishra}
\seealso{
\code{\link{calWt}},
\code{\link{RAWgrid}},
\code{\link{nb}},
\code{\link{cvWtTuning}}
}
\examples{
### Load data ##
\dontrun{
data(fakeData)

### Get grid of tuning parameters  ###
grid <- RAWgrid(r = 0.3,rl = -Inf,ru = Inf,p = fakeData$p,y = fakeData$y,
                cvParm = "lambda",rl.raw = 0.25,ru.raw = 0.35)

### Implement repeated k-fold cross validation
repCV <- cvRepWtTuning(y = fakeData$y,p = fakeData$p,rl = -Inf,ru = Inf,r = 0.3,
                       kFold = 5,cvRep = 25,cvParm = "lambda",tuneSeq = grid,stdErrRule = TRUE)

## cross-validation results
repCV$avgCV.res

## cross-validation selected lambda, RAW, and sNV
cv.lambda <- repCV$cv.lambda
cv.RAW <- repCV$cv.RAW
cv.RAW <- repCV$cv.sNB
}
}
