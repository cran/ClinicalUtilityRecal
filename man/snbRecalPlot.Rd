\name{snbRecalPlot}
\alias{snbRecalPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of Potential sNB Under Recalibration
}
\description{
Function for producing a graphical device to help assess the potential for recalibration to improve the clinical utility of a risk model.

This function plots the sNB for a given risk score (p) and cost benefit ratio (R/1-R) against the threshold used to perscribe intervention. The peak of the curve estimates the maxmium possible sNB that can be obtained via recalibration. The estimated sNB of the original risk score and the recalibrated risk score under standard logistic recalibration are also plotted on the curves on the curves. The dotted lined  "stdErrThresh" controls. The plot includes a dotted horizontal line one standard error (or any number of standard errors, controlled by stdErrThresh) below the maximum, to help guage much room for improvement in sNB there is for original or standard recalibrated risk score from maximum possible sNB.
}
\usage{
snbRecalPlot(p,p.std,y,r,stdErrThresh=1,ylim=NULL,
              titlePlot = "Potential sNB Under Recalibration",risk.model.std=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{Vector of risk score values}
  \item{p.std}{Vector of risk score values after standard logistic recalibration}
  \item{y}{Vector of binary outcomes, with 1 indicating event (case) and 0 indicating no event (controls)}
  \item{r}{Clinically relevant risk threshold}
  \item{stdErrThresh}{Indicates how many standard errors line drawn below the maximum of the sNB curve should be}
  \item{ylim}{Y axis limits}
  \item{titlePlot}{Title for plot}
  \item{risk.model.std}{Plot standard error bars for sNB of plotted risk models}

}

\value{Prints plot of potential sNB of risk score p under recalibration}


\references{Mishra, A. (2019). Methods for Risk Markers that Incorporate Clinical Utility (Doctoral dissertation). (Available Upon Request)}

\author{
Anu Mishra
}


\examples{
data("fakeData")

### get standard recalibrated risk score
stdRecal.res <- stdRecal(y = fakeData$y,p = fakeData$p)
p.std <- stdRecal.res$p.std


## make plot
snbRecalPlot(p = fakeData$p,p.std = p.std,y = fakeData$y,r = 0.3)

## both original and std logistic recalibrated risk score are below 1 std err from maximum
## indicating that alternative recalibration methods could improve
}
