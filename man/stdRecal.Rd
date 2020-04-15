\name{stdRecal}
\alias{stdRecal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard Logistic Recalibration
}
\description{
Produces recalibration intercept, slope and corresponding recalibrated risk scores using standard logistic rexalibration. Recalibration slope and intercept using logistic recalibration method, develop by Cox (1958), then used to calculated the recalibrated risk score
}
\usage{
stdRecal(y,p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{
  Vector of binary outcomes, with 1 indicating event (case) and 0 indicating no event (controls)}
   \item{p}{
  Vector of risk score values}
}
\value{
\item{stdRisk }{Vector of recalibrated risks under standard logistic recalibration }
\item{alpha }{Two-element vector containing logistic recalibration intercept and slope}
}
\references{
Cox, D. R. (1958). Two further applications of a model for binary regression. Biometrika, 45(3/4), 562-565.

Harrell Jr, F. E. (2015). Regression modeling strategies: with applications to linear models, logistic and ordinal regression, and survival analysis. Springer.
}

\author{
Anu Mishra
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
#load data
data(fakeData)

res <- stdRecal(p=fakeData$p, y =fakeData$y)
p.std <- res$stdRecal
coef.dt <- res$alpha
}
