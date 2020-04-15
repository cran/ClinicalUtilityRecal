\name{constRecal}
\alias{constRecal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constrained Logistic Recalibration
}
\description{
When recalibrating a risk model, where the intended purpose of the risk model is to prescribe an intervention to those deemed at high risk, it is desirable to have good calibration at the clinically relevant threshold used to define high risk (i.e. threshold used to identify who recieves treatment). This can be achieved by maximizing the clinical utilty of a risk model, which should in turn produce good calibration near the clinical relevant risk threshold. This function realibrates risk scores (predicting binary outcome) by estimating the recalibration intercept and slope by maximizing the logistic likelihood over a restricted parameter space (Mishra et al. [2020]).

The restricted space only includes recalibration parameters \eqn{alpha_0} and \eqn{alpha_1} that produce a recalibrated risk model with high sNB. The restricted parameter space is defined all \eqn{alpha} parameters that result in a recalibrated risk score within one-standard-error of the maximum possible \eqn{sNB}. See Mishra et al (2020), for full details}

\usage{
constRecal(y,p,r,int=NULL,alphaLB=c(-10,0),
            alphaUB=c(10,10),ftol=1e-8,xtol=1e-4,maxeval=1e6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{
  Vector of binary outcomes, with 1 indicating event (case) and 0 indicating no event (controls)}
   \item{p}{
  Vector of risk score values}
      \item{r}{
  Clinically relevant risk threshold}
  \item{int}{Two-dimensional vector of initial recalibration parameter (\eqn{alpha}) values. If not specified the recalibration parameters estimated under standard logistic recalibration will be used}
 \item{alphaLB}{Lower bound of box-contrained search space}
  \item{alphaUB}{Upper bound of box-contrained search space}
  \item{ftol}{Controls tolerance of optimization procedure with respect to changes in sNB}
  \item{xtol}{Controls tolerance of optimization procedure with respect to changes in alpha}
  \item{maxeval}{Maximum number of interations performed during optimization procedure}
}
\details{
To solve this optimization problem the DIRECT optimization method is implemented via the NLOPTR package. See Jones et al (1993) and Ypma et al (2014) for full description of optimization method and implementation details. Note this is not a convex optimization problem, so a global optimizer is used.
}
\value{
\item{alpha}{Recalibration parameters obtained from Constrained Logistic Recalibration}
\item{y}{Vector of binary outcomes, with 1 indicating event (case) and 0 indicating no event (controls)}
\item{p.const}{Vector of constrained logistic reclaibrated risk scores}
\item{snbLB}{Lower bound used to define constraint region. The solution is constrained to the region where recalibration parameters result in risk model with \eqn{sNB} greater or equal to this bound. This is one-standard error lower than the maximum \eqn{sNB} }
\item{optimRes}{Full output from the NLOPTR optimization routine}
}
\references{Mishra, A. (2019). Methods for Risk Markers that Incorporate Clinical Utility (Doctoral dissertation). (Available Upon Request)

Ypma, J., Borchers, H. W., Eddelbuettel, D., & Ypma, M. J. (2020). Package ‘nloptr’.

Ypma, J. (2014). Introduction to nloptr: an R interface to NLopt. Tech. rep.

D. R. Jones, C. D. Perttunen, and B. E. Stuckmann, "Lipschitzian optimization without the lipschitz constant," J. Optimization Theory and Applications, vol. 79, p. 157 (1993).

}
\author{
Anu Mishra
}

\examples{
\dontrun{
### Load data ##
data(fakeData)

## Implementing standard logistic recalibration
stdRecal.res <- stdRecal(y = fakeData$y,p = fakeData$p)
stdRecal.res$alpha #standard recalibration parameters
p.std <- stdRecal.res$p.std

## Look at potential sNB under recalibration plot
snbRecalPlot(p = fakeData$p,p.std = p.std,y = fakeData$y,r = 0.3)

## Implementing constrained logistic recalibration
constRecal.res <- constRecal(y = fakeData$y,p = fakeData$p,r = 0.3)
constRecal.res$alpha #constrained logistic recalibration parameters
p.recal <- constRecal.res$p.const

## comparing standardized net benefit of the two
nb(y = fakeData$y,p = fakeData$p,r = 0.3)$snb #original
nb(y = stdRecal.res$y,p = stdRecal.res$p.std,r = 0.3)$snb #std recal
nb(y = constRecal.res$y,p = constRecal.res$p.const,r = 0.3)$snb #weighted

## Generate calibration plots

### Calibration curve of only original, standard and weighted recalibrated risk score
calCurvPlot(y=fakeData$y,p=fakeData$p,p.std=p.std,p.recal=p.recal,
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
