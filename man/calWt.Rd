\name{calWt}
\alias{calWt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calibration Weights
}
\description{
Calculates observation weights used for weighted calibration method using LOESS smoother. Observations with predicted risks outside clinically relevant interval [Rl,Ru] are downweighted.
}
\usage{
calWt(rl,ru,p,y,r,lambda,delta,returnSmoothedEvent=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rl}{
  Lower bound of clinically relevant region}
    \item{ru}{
  Upper bound of clinically relevant region}
    \item{p}{
  Vector of risk score values}
    \item{y}{
  Vector of binary outcomes, with 1 indicating event (cases) and 0 indicating no event (controls)}
    \item{r}{
  Clinically relevant risk threshold}
    \item{lambda}{
  Tuning parameter for weights inside relevant region}
  \item{delta}{
  Weight assigned to observations outside relevant region
  }
  \item{returnSmoothedEvent}{
  If TRUE returns smoothed observed event rate used for calculationg weights
  }

}
\details{
Computes observation weights for weighted recalibration method. Observations with risk scores near a pre-defined clinically relevant risk threshold \eqn{r} are given weights near 1, while observations far from the risk threshold are down-weighted.

For observations with risk scores within a pre-defined clinically relevant region [R_l, R_u], the weighting function follows an exponential decay form with observations recieving smaller weights as they move farther from the clincally relevant risk threhsold. Distance is measured by the squared difference between the risk threshold and \eqn{o(p_i)}, where \eqn{o(p_i)} isa smoothed observed event rate,obtained via LOESS regression of \eqn{y} on the risk scores \eqn{p_i}.

The amount of decay is set by the tuning parameter \eqn{\lambda}, with larger lambda indicating little down-weighting applied. For observations outside the clinically relevant risk interval, a weight of \eqn{delta} is assigned. An indicator type weight can be choosen by select large \eqn{\lambda} (e.g \eqn{\lambda=10}). Clinically relevant region can be [R_l, R_u] = [\eqn{-\infty},\eqn{\infty}], and in these cases it is not necessary to specify \eqn{delta}. See Mishra et al (2020) for functional form of weights and more details.}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{wt}{Vector of weights, between 0 and 1, to be used for weighted recalibration}
 \item{o}{Vector of smoothed observed event rates}
}

\references{Mishra, A. (2019). Methods for Risk Markers that Incorporate Clinical Utility (Doctoral dissertation). (Available Upon Request)}




\author{
Anu Mishra
}

\seealso{
\code{\link{wtRecal}}
}

\examples{
#load data
data(fakeData)

#exponential decay type weight
expWt <- calWt(rl=-Inf,ru=Inf,p=fakeData$p,y=fakeData$y,r=0.2,
                lambda=0.05,delta=0.2,returnSmoothedEvent = TRUE)
plot(expWt$o,expWt$wt,xlab="Smoothed Observed Event Rate",
     ylab="Calibration Weight",ylim=c(0,1))


#indicator type weight
indWt <- calWt(rl=0.1,ru=0.3,p=fakeData$p,y=fakeData$y,r=0.2,
                lambda=0.05,delta=0.2,returnSmoothedEvent = TRUE)
plot(indWt$o,indWt$wt,xlab="Smoothed Observed Event Rate",
    ylab="Calibration Weight",ylim=c(0,1))
}

