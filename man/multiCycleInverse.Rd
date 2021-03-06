% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PsiFunctions.R
\name{multiCycleInverse}
\alias{multiCycleInverse}
\title{multiCycleInverse}
\usage{
multiCycleInverse(warpedStimuli, references = c(0))
}
\arguments{
\item{warpedStimuli}{a vector of stimuli, between -inf and inf, or a list of vectors of stimuli}

\item{references}{A vector of values.  This should include -Inf and Inf, but these are implicitly added
if they are omitted}
}
\value{
A vector containing warped stimuli, or a list of vectors
}
\description{
recombines a mapping, of any sort (bounded, semi-bounded, or unbounded) into single region,
by eliminating 'references' inside the space. Functionally, multiCycleInvers simply undoes a vector of stimuli into
}
\details{
For instance, a spatial region |----------------|
might be divided by a line of vertical symmetry into two regions
 |--------||--------|
 And each be turned int an unbounded Prelec region
 <--------><-------->
 with a command like -10:10 %>% multicycle(0) %>% psiPrelec
 multiCycleInverse would then recombine these values into a single range.

KEY LIMITATION: Right now, multiCycleInverse can only 'handle' one layer of recursion. That is, you can't make calls like
 -10:10 %>% multiCycle %>% multiCycle %> multiCycleInverse %>% multiCycleInverse
 . Sorry.
}
\examples{
multiCycle(-99:100, c(-100, 0, 100))
plot(-99:100, unlist(multiCycle(-99:100, c(-100, -50, 0, 50, 100))))
(-99:100/100) \%>\% multiCycle(-1, 0, 1) \%>\% 
    psiLogOdds() \%>\% vanillaBayes() \%>\% 
    psiLogOddsInverse() \%>\% 
    multiCycleInverse(-1, 0, 1) # Implements Landy et al's model 
    # of one-dimensional spatial memory, with fixed boundaries
}
\seealso{
psiIdentity, psiLogOdds, psiPrelec, psiLog, psiLinearInverse
}
\keyword{perceptual}
\keyword{psi}
\keyword{tranformations}
