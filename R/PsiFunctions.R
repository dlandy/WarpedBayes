# Various Psi functions, along with their inverses, with the following naming convention
# psiName <- function(stimuli, parameters with defaults)
# psiNameInverse <- function(unboundedSpace, parameters matching psiName)

#' psiLinear
#' 
#' Creates mappings from -inf:+inf -> -inf:+inf, with a linear transformation mechanism
#' @param stimuli a vector of stimuli, between -inf and inf
#' @return a vector containing warped stimuli
#' @keywords psi, perceptual tranformations
#' @seealso psiIdentity, psiLogOdds, psiPrelec, psiLog
#' @export
#' @examples
#' psiLinear(c(0.1, 0.2, 0.3))
#' psiLinear(-10:10)
#' psiLinear(-10:10, shift=10, scaling=2)
psiLinear <- function(stimuli, shift=0, scaling=1){
  stimuli*scaling + shift
}



#' psiIdentity
#' 
#' This function is a generic transformation function that does no conversion at all. 
#' It's really just here for testing things out.
#' It is equivalent to calling psiLinear with default parameters
#' @param stimuli a vector of stimuli, between -inf and inf
#' @return a vector containing "warped" stimuli identical to the original stimuli
#' @keywords psi, perceptual tranformations
#' @seealso psiIdentity, psiLogOdds, psiPrelec, psiLog
#' @export
#' @examples
#' psiIdentity(-10:10)
psiIdentity <- function(stimuli){
  stimuli
}
