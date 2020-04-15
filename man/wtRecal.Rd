\name{wtRecal}
\alias{wtRecal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Logistic Recalibration
}
\description{
When recalibrating a risk model, where the intended purpose of the risk model is to prescribe an intervention to those deemed at high risk, it is desirable to have good calibration at the clinically relevant threshold used to define high risk (i.e. threshold used to identify who recieves treatment). This function realibrates risk scores (predicting binary outcome) using the weighted logistic recalibration method (Mishra et al. [2020]).

Under this method, a recalibration intercept and slope are estimated via weighted logistic regression. Weights are constructed such that observations further from the clinically relevant risk threshold are down-weighted, meaning observations closer to the clinically relevant risk threshold have higher contribution to the risk threshold. The resulting estimated weighted recalibration slope and intercept are used to scale and shift the existing risk score, producing better calibrated risk scores near the risk threhold and potentially increaseing the net benefit of the risk score.
}
\usage{
wtRecal(y,p,r,rl,ru,lambda,delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{
  Vector of binary outcomes, with 1 indicating event (case) and 0 indicating no event (controls)}
   \item{p}{
  Vector of risk score values}
      \item{r}{
  Clinically relevant risk threshold}
    \item{rl}{
  Lower bound of clinically relevant region}
    \item{ru}{
  Upper bound of clinically relevant region}
    \item{lambda}{
  Tuning parameter for weights inside relevant region}
  \item{delta}{
  Weight assigned to observations outside relevant region
  }
}

\value{
\item{y}{Vector of binary outcomes, with 1 indicating event (case) and 0 indicating no event (controls)}
\item{p.wt}{Vector of weighted reclaibrated risk scores}
\item{alpha.wt}{Weighted recalibration slope and intercept}
\item{wt}{Calibration weights}
\item{wt.conv}{Indicator of convergence of the weighted logistic regression model. 1 indicates model convergence, 0 indicates model did not converge}
}
\references{Mishra, A. (2019). Methods for Risk Markers that Incorporate Clinical Utility (Doctoral dissertation). (Available Upon Request)}


\author{
Anu Mishra
}
\seealso{
\code{\link{calWt}}
}
\examples{
\dontrun{
### Load data ##
data(fakeData)

## Implement standard logistic recalibration
stdRecal.res <- stdRecal(y = fakeData$y,p = fakeData$p)
stdRecal.res$alpha #standard recalibration parameters
p.std <- stdRecal.res$p.std

## Look at potential sNB under recalibration plot
snbRecalPlot(p = fakeData$p,p.std = p.std,y = fakeData$y,r = 0.3)
## both original and std logistic recalibrated risk score are below 1 std err from maximum
## indicating that alternative recalibration methods could improve

### Get grid of tuning parameters  ###
## In this example keeping delta fixed, and selecting lambda

grid <- RAWgrid(r = 0.3,rl = -Inf,ru = Inf,p = fakeData$p,y = fakeData$y,
                cvParm = "lambda",rl.raw = 0.25,ru.raw = 0.35)


### Select tuning parameter lamba using 5-fold cross-validation repeated 25 times
### with one standard error rule implemented

repCV <- cvRepWtTuning(y = fakeData$y,p = fakeData$p,rl = -Inf,ru = Inf,r = 0.3,
                       kFold = 5,cvRep = 25,cvParm = "lambda",tuneSeq = grid,stdErrRule = TRUE)

## Implement weighted logistic recalibration
wtRecal.res <- wtRecal(y = fakeData$y,p = fakeData$p,r = 0.3,rl = -Inf,ru = Inf,
                       lambda = repCV$cv.lambda,delta=1)
# note that delta here is set to 1 as a place holder, but not used since the clinically
# relevant region is [-Inf, Inf]
wtRecal.res$alpha.wt #weighted recalibration parameters
p.recal <- wtRecal.res$p.wt

## comparing standardized net benefit of the two
nb(y = fakeData$y,p = fakeData$p,r = 0.3)$snb #original
nb(y = stdRecal.res$y,p = stdRecal.res$p.std,r = 0.3)$snb #std recal
nb(y = wtRecal.res$y,p = wtRecal.res$p.wt,r = 0.3)$snb #weighted

### Calibration curve of only original, standard and weighted recalibrated risk score
calCurvPlot(y = fakeData$y,p = fakeData$p,p.std=p.std,p.recal=p.recal,
            stdPlot=TRUE, recalPlot=TRUE,
            xlim=c(0,1),ylim=c(0,1),
            label="Original Risk Score",
            label2 = "Standard Recalibrated Risk Score",
            label3 = "Weighted/Constrained Recalibrated Risk Score",
            legendLab = c("Orig.", "Std.", "Wt."),
            mainTitle="Calibration of Risk Score",
            hist=TRUE,ylimHist = c(0,0.5),
            r=0.3,rl = -Inf, ru = Inf)
}

}
