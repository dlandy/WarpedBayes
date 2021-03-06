% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PsiFunctions.R
\name{psiLogOdds}
\alias{psiLogOdds}
\title{psiLogOdds}
\usage{
psiLogOdds(stimuli, smallValue = 10^-30)
}
\arguments{
\item{stimuli}{a vector of positive stimuli, between 0 and inf}

\item{smallValue}{A value to set nominal 0's to, to avoid errors in plotting simulated data, or aberrant responses}
}
\value{
a vector containing "warped" stimuli
}
\description{
This function is a transformation function that takes the log odds
By design, it is robust to the inclusion of the occasional 0,
which it maps to a very small value, which can optionally be specified
It's appropriate for mapping (0, 1):-> (-inf, inf).
}
\examples{
psiLogOdds(1:100/100)
(0:1000/1000) \%>\% psiLogOdds() \%>\%
     vanillaBayes() \%>\% 
     psiLogOddsInverse()  # Implements Gonzales & Wu, 1996
}
\seealso{
psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
}
\keyword{perceptual}
\keyword{psi}
\keyword{tranformations}
