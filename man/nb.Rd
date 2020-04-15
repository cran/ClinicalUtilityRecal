\name{nb}
\alias{nb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Net Benefit
}
\description{
Estimates the net benefit and standardized net benefit of a risk score given a risk threshold R. Additionally estimates the standardized net benefit of a "treat-all" (everyone recieves intervention) treatment rule.
}
\usage{
nb(y, p, r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{
  Vector of binary outcomes, with 1 indicating event (case) and 0 indicating no event (controls)}
   \item{p}{
  Vector of risk score values}
      \item{r}{
  Clinically relevant risk threshold}
}
\value{
\item{nb }{Net benefit of risk score \eqn{p}}
\item{snb }{Standardized net benefit of risk score \eqn{p}. Net benefit of \eqn{p} dived by the prevalence}
\item{snb.all}{Standardized net benefit of a decision rule where everyone recieves treatment }
}

\references{
Pauker, S. G., & Kassirer, J. P. (1980). The threshold approach to clinical decision making. New England Journal of Medicine, 302(20), 1109-1117.

Mishra, A. (2019). Methods for Risk Markers that Incorporate Clinical Utility (Doctoral dissertation). (Available Upon Request)}

\author{
Anu Mishra }
