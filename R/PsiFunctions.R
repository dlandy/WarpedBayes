# Various Psi functions, along with their inverses, with the following naming convention
# psiName <- function(stimuli, parameters with defaults)
# psiNameInverse <- function(unboundedSpace, parameters matching psiName)

#' psiLinear
#' 
#' Creates mappings from -inf:+inf -> -inf:+inf, with a linear transformation mechanism
#' @param stimuli a vector of stimuli, between -inf and inf, or a list of vectors of stimuli
#' @param shift A scalar by which to (positively) offset the scaled values, or a list or vector of scalars
#' @param scaling A scalar multiplier to the (unshifted) values, or a list or vector of scalars
#' @return A vector containing warped stimuli, or a list of vectors
#' @keywords psi perceptual tranformations
#' @seealso psiIdentity, psiLogOdds, psiPrelec, psiLog, psiLinearInverse
#' @export
#' @examples
#' psiLinear(c(0.1, 0.2, 0.3))
#' psiLinear(-10:10)
#' psiLinear(-10:10, shift=10, scaling=2)
#' -10:10 %>% psiLinear(shift=2, scaling=0.5)
psiLinear <- function (x, ...) {
  UseMethod("psiLinear", x)
}


psiLinear.double <- function(stimuli, shift=0, scaling=1){
    stimuli*scaling + shift
}

psiLinear.list <- function(stimuli, shift=0, scaling=1){
    mapply("+", mapply("*", stimuli, scaling, SIMPLIFY = FALSE), shift, SIMPLIFY=FALSE)
}
  

#' psiIdentity
#' 
#' This function is a generic transformation function that does no conversion at all. 
#' It's really just here for testing things out.
#' It is equivalent to calling psiLinear with default parameters
#' @param stimuli a vector of stimuli, between -inf and inf
#' @return a vector containing "warped" stimuli identical to the original stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiIdentity(-10:10)
psiIdentity <- function(stimuli){
  stimuli
}


#' psiLog
#' 
#' This function is a transformation function that takes the log...
#' It's literally just log(stimuli).
#' It's appropriate for mapping (0, +inf):-> (-inf, inf).
#' By design, it is robust to the inclusion of the occasional 0, which it maps to a very small value.
#' @param stimuli a vector of positive stimuli, between 0 and inf
#' @param smallValue A value to set nominal 0's to, to avoid errors in plotting simulated data, or aberrant responses
#' @return a vector containing "warped" stimuli identical to the original stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiLog(1:100)
#' 1:100 %>% psiLog()
#' 1:1000 %>% psiLog() %>% vanillaBayes() %>% psiLogInverse()  # Implements Stevens Power Law

psiLog <- function(stimuli, smallValue=10^-5){
  stimuli[stimuli==0] <- smallValue
  log(stimuli)
}

#' psiLogOdds
#' 
#' This function is a transformation function that takes the log odds
#' By design, it is robust to the inclusion of the occasional 0,
#' which it maps to a very small value, which can optionally be specified
#' It's appropriate for mapping (0, 1):-> (-inf, inf).
#' @param stimuli a vector of positive stimuli, between 0 and inf
#' @param smallValue A value to set nominal 0's to, to avoid errors in plotting simulated data, or aberrant responses
#' @return a vector containing "warped" stimuli identical to the original stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiLogOdds(1:100/100)
#' 0:1000/1000 %>% psiLogOdds() %>% vanillaBayes() %>% psiLogOddsInverse()  # Implements Gonzales & Wu, 1996
psiLogOdds <- function(stimuli, smallValue=10^-5){
  stimuli[stimuli==1] <- 1-smallValue
  stimuli[stimuli==0] <- smallValue
  d <- log(stimuli/(1-stimuli))
  d
  }



#INVERSE Functions



#' psiIdentityInverse
#' 
#' This function is a generic transformation function that does no conversion at all. 
#' It is the inverse of psiIdentity.
#' It's really just here for testing things out.
#' It is equivalent to calling psiLinearInverse with default parameters
#' @param warpedStimuli a vector of stimuli, between -inf and inf
#' @return a vector containing "unwarped" stimuli identical to the warped stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentity
#' @export
#' @examples
#' psiIdentity(-10:10)
psiIdentity <- function(warpedStimuli){
  warpedStimuli
}



#' psiLinearInverse
#' 
#' Creates mappings from -inf:+inf -> -inf:+inf, with a linear transformation mechanism
#' @param warpedStimuli a vector of stimuli, between -inf and inf
#' @return A vector containing unwarped stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiIdentity, psiLogOdds, psiPrelec, psiLog, psiLinear
#' @export
#' @examples
#' psiLinear(c(0.1, 0.2, 0.3))
#' psiLinear(-10:10)
#' psiLinear(-10:10, shift=10, scaling=2)
#' -10:10 %>% psiLinear(shift=2, scaling=0.5)
psiLinearInverse <- function(warpedStimuli, shift=0, scaling=1){
  (warpedStimuli-shift)/scaling 
}



#' psiLogInverse
#' 
#' This function is a transformation function that inverts the log...
#' It's literally just exp(stimuli).
#' It's appropriate for mapping  (-inf, inf):-> (0, +inf)
#' @param stimuli a vector of positive stimuli, between 0 and inf
#' @return a vector containing "warped" stimuli identical to the original stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiLogInverse(psiLog(1:100)) # returns 1:100
#' psiLog(1:100) %>% psiLogInverse()
#' 1:1000 %>% psiLog() %>% vanillaBayes() %>% psiLogInverse()  # Implements Stevens Power Law

psiLogInverse <- function(warpedStimuli, smallValue=10^-5){
  exp(warpedStimuli)
}

#' psiLogOddsInverse
#' 
#' This function is a transformation function that inverts the log odds

#' It's appropriate for mapping (-inf, inf) :->  (0, 1)
#' As in all functions in this package, the parameters are set so that applying the same parameters to the main function 
#' and the inverse yields an identity.
#' @param stimuli a vector of positive stimuli, between 0 and inf
#' @param smallValue A value to set nominal 0's to, to avoid errors in plotting simulated data, or aberrant responses. Accepted here to give the function and its inverse the same parameter set
#' @return a vector containing "warped" stimuli identical to the original stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiLogOdds(1:100/100)
#' 0:1000/1000 %>% psiLogOdds() %>% vanillaBayes() %>% psiLogOddsInverse()  # Implements Gonzales & Wu, 1996
psiLogOddsInverse <- function(warpedStimuli, smallValue=10^-5){
  e <- exp(warpedStimuli)
  d <- e/(1+e)
  d
}



# MULTICYCLING

# Multicycling is an advanced art.  Basically, we 'multicycle' by a pair of functions and inverses 
# that carve spaces into regions that are already bounded at lynchpin points.  These can include 'inf' and '-inf'



