% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VanillaBayes.R
\name{vanillaBayes}
\alias{vanillaBayes}
\title{vanillaBayes}
\usage{
vanillaBayes(stimuli, kappa = 0, tauStimuli = 1, tauCategory = 1,
  responses = NULL, mode = "prediction", responseGrid = c(0))
}
\arguments{
\item{stimuli}{a vector of stimuli, between -inf and inf}

\item{kappa}{The location of the category}

\item{tauStimuli}{The precision of the stimulus traces: may be a single number or a vector}

\item{tauCategory}{The precision of the category distribution}

\item{responses}{an optional vector of responses. Should only be given if mode is "subjectiveLogLikelihood"}

\item{mode}{What aspect should the function calculate? Legel choices include "prediction", "simulation", and "subjectiveLogLikelihood"}
}
\value{
A vector containing mean predicted stimulus locations, or the log likelihood of the responses given the model
}
\description{
Performs vanilla bayes (single category) on a set of stimuli.
}
\details{
This function assumes that the data are in a metric 
space (-inf, inf), with a single normally distributed generating category (with mean kappa and precision tauCategory). It further assumes a set of stimuli, which are
normal distributions with means at the value of stimuli, and precision tauStimuli.  It returns either the mean expected 
location of response to the stimuli (given the parameters), or if a set of responses is given, the log likelihood of the 
responses given the model.

If the kappa, tauStimuli, and tauCategory items are all more than length 1,
 and  are length 2 less than the number of bins, then we pad them by negative and positive infinity.
}
\examples{
(0:1000/1000) \%>\% vanillaBayes(kappa=5) 
    # The  Bayesian normal-normal model typical to many analyses
(0:1000/1000) \%>\% psiLogOdds() \%>\% 
    vanillaBayes(kappa=5) \%>\% 
    psiLogOddsInverse()  #Gonzales & Wu, 1996
1:1000 \%>\% psiLog() \%>\% vanillaBayes() \%>\% psiLogInverse()  # Stevens Power Law
plot(-99:100/100, (-99:100/100) \%>\% multiCycle(references= c(-10, 0, 10)) \%>\% 
    psiLogOdds() \%>\% vanillaBayes(kappa=c(-1, 1), tauStimuli=10) \%>\% 
    psiLogOddsInverse() \%>\% 
    multiCycleInverse(references=c(-10, 0, 10))-(-99:100/100), 
        ylab="bias", xlab="stimulus");abline(0,0)
}
\seealso{
vanillaBayes
}
\keyword{bayesianInference}
