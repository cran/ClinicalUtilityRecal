\name{calCurvPlot}
\alias{calCurvPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calibration Curve Plot
}
\description{
Flexible function for plotting the calibration curve for a risk score and different recalibrations of a risk score. Histograms of risk score and recalibrated risk scores can be plotted using argument "hist==TRUE". Dotted guidelines for the risk threshold R are plotted and clinically releavant risk region [R_l, R_u] are plotted to help gauge calibration in clinically relevants.
}
\usage{
calCurvPlot(y,p,p.std=NULL,p.recal=NULL,
            stdPlot=FALSE, recalPlot=FALSE,
            xlim=c(0,1),ylim=c(0,1),
            label="Original Risk Score",
            label2 = "Standard Recalibrated Risk Score",
            label3 = "Weighted/Constrained Recalibrated Risk Score",
            legendLab = c("Orig.", "Std.", "Wt."),
            mainTitle="Calibration of Risk Score",
            hist=TRUE,ylimHist = c(0,0.5),
            r,rl = -Inf, ru = Inf)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Vector of binary outcomes, with 1 indicating event (case) and 0 indicating no event                 (controls)}
  \item{p}{Vector of risk score values}
  \item{p.std}{Vector of risk score values after standard logistic recalibration}
  \item{p.recal}{Vector of risk score values after weighted/contratined logistic recalibration}
  \item{stdPlot}{If TRUE plot calibration curve for standard logistic recalibrated risk score }
  \item{recalPlot}{If TRUE plot calibration curve for weighted/contratined recalibrated risk score }
  \item{xlim}{Limits for x-axis}
  \item{ylim}{Limits for y-axis}
  \item{label}{Label for x-axis corresponding to p vector}
  \item{label2}{Label for x-axis corresponding to p.std vector}
  \item{label3}{Label for x-axis corresponding to p.recal vector}
  \item{legendLab}{Label for legend}
  \item{mainTitle}{Main title for plot}
  \item{hist}{If true plot distribution of risk scores along with calibration curve}
  \item{ylimHist}{Limits for y-axis of histogram}
  \item{r}{Clinically relevant risk threshold used indicated with dotted line}
  \item{rl}{Lower bound of clinically relevant region indicated with dotted-dashed line}
  \item{ru}{Upper bound of clinically relevant region indicated with dotted-dashed line}
}

\value{Prints calibration plot of risk model p, and optional histogram of risk score}

\references{Mishra, A. (2019). Methods for Risk Markers that Incorporate Clinical Utility (Doctoral dissertation). (Available Upon Request)}

\author{
Anu Mishra
}

\examples{
\dontrun{
data("fakeData")

### get standard recalibrated risk score
stdRecal.res <- stdRecal(y = fakeData$y,p = fakeData$p)
p.std <- stdRecal.res$p.std

### Select tuning parameter lamba using 5-fold cross-validation repeated 25 times
grid <- RAWgrid(r = 0.3,rl = -Inf,ru = Inf,p = fakeData$p,y = fakeData$y,
                cvParm = "lambda",rl.raw = 0.25,ru.raw = 0.35)
repCV <- cvRepWtTuning(y = fakeData$y,p = fakeData$p,rl = -Inf,ru = Inf,r = 0.3,
                       kFold = 5,cvRep = 25,cvParm = "lambda",tuneSeq = grid,stdErrRule = TRUE)

## Implement weighted logistic recalibration
wtRecal.res <- wtRecal(y = fakeData$y,p = fakeData$p,r = 0.3,rl = -Inf,ru = Inf,
                       lambda = repCV$cv.lambda,delta=1)
p.recal <- wtRecal.res$p.wt

### Calibration curve of only original risk score with histogram
calCurvPlot(y=fakeData$y,p=fakeData$p,p.std=NULL,p.recal=NULL,
            stdPlot=FALSE, recalPlot=FALSE,
            xlim=c(0,1),ylim=c(0,1),
            label="Original Risk Score",
            label2 = "Standard Recalibrated Risk Score",
            label3 = "Weighted/Constrained Recalibrated Risk Score",
            legendLab = c("Orig.", "Std.", "Wt."),
            mainTitle="Calibration of Risk Score",
            hist=TRUE,ylimHist = c(0,0.5),
            r=0.3,rl = -Inf, ru = Inf)

### Calibration curve of only original, standard and weighted recalibrated risk score
calCurvPlot(fakeData$y,p=fakeData$p,p.std=p.std,p.recal=p.recal,
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
