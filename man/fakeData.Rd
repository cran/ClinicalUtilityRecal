\name{fakeData}
\alias{fakeData}
\docType{data}
\title{
Dataset for Recalibration Purposes
}
\description{
Simulated dataset that can be used for recalibration purposes.
}
\usage{data("fakeData")}
\format{
  A data frame with 1000 observations on the following 2 variables.
  \describe{
    \item{\code{p}}{a numeric vector or risk scores bounded between 0 and 1}
    \item{\code{y}}{a numeric vector of indicator of events (cases) and non-events (control)}
  }
}
\details{
See Mishra et al (2020) for description of simuation settings to generate data (simulation example 3).
}

\references{Mishra, A. (2019). Methods for Risk Markers that Incorporate Clinical Utility (Doctoral dissertation). (Available Upon Request)}


\examples{
data(fakeData)
}
\keyword{datasets}
