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
psiLinear <- function (stimuli, shift=0, scaling=1) {
  UseMethod("psiLinear",stimuli)
}

#' @export
psiLinear.list <- function(stimuli, shift=0, scaling=1){
    mapply(psiLinear, stimuli, shift, scaling, SIMPLIFY=FALSE)
}

#' @export
psiLinear.subjectiveLogLikelihood <- function(stimuli, shift=0, scaling=1){  stimuli}

#' @export
psiLinear.default <- function(stimuli, shift=0, scaling=1){stimuli}

#' @export
psiLinear.numeric <- function(stimuli, shift=0, scaling=1){
  stimuli*scaling + shift
}  

#' psiIdentity
#' 
#' This function is a generic transformation function that does no conversion at all. 
#' It's really just here for testing things out.
#' It is equivalent to calling psiLinear with default parameters
#' @param stimuli a vector of stimuli, between -inf and inf
#' @param smallValue A value to set nominal 0's to, to avoid errors in plotting simulated data, or aberrant responses
#' @return a vector containing "warped" stimuli identical to the original stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiIdentity(-10:10)
psiIdentity <- function (stimuli, smallValue=10^-30) {
  UseMethod("psiIdentity", stimuli)
}

#' @export
psiIdentity.list <- function(stimuli, smallValue=10^-30){
  mapply(psiIdentity, stimuli, SIMPLIFY=FALSE)
}

#' @export
psiIdentity.default <- function(stimuli, smallValue=10^-30){stimuli}


#' @export
psiIdentity.subjectiveLogLikelihood <- function(stimuli, smallValue=10^-30){stimuli}


#' @export
psiIdentity.numeric <- function(stimuli, smallValue=10^-30){
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
#' @return a vector containing "warped" stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiLog(1:100)
#' 1:100 %>% psiLog()
#' 1:1000 %>% psiLog() %>% 
#'     vanillaBayes() %>% 
#'     psiLogInverse()  # Implements Stevens Power Law
psiLog <- function (stimuli, smallValue=10^-30) {
  UseMethod("psiLog", stimuli)
}

#' @export
psiLog.list <- function(stimuli, smallValue=10^-30){
  mapply(psiLog, stimuli, smallValue, SIMPLIFY=FALSE)
}

#' @export
psiLog.subjectiveLogLikelihood <- function(stimuli, smallValue){stimuli}


#' @export
psiLog.default <- function(stimuli, smallValue){stimuli}



#' @export
psiLog.numeric <- function(stimuli, smallValue=10^-30){
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
#' @return a vector containing "warped" stimuli 
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiLogOdds(1:100/100)
#' (0:1000/1000) %>% psiLogOdds() %>%
#'      vanillaBayes() %>% 
#'      psiLogOddsInverse()  # Implements Gonzales & Wu, 1996
psiLogOdds <- function (stimuli, smallValue=10^-30) {
  UseMethod("psiLogOdds", stimuli)
}

#' @export
psiLogOdds.list <- function(stimuli, smallValue=10^-30){
  mapply(psiLogOdds, stimuli, smallValue, SIMPLIFY=FALSE)
}

#' @export
psiLogOdds.default <- function(stimuli, smallValue){stimuli}



#' @export
psiLogOdds.subjectiveLogLikelihood <- function(stimuli, smallValue){stimuli}


#' @export
psiLogOdds.numeric <- function(stimuli, smallValue=10^-30){
  stimuli[stimuli>= 1-smallValue] <- 1-smallValue
  stimuli[stimuli <= smallValue] <- smallValue
  d <- log(stimuli/(1-stimuli))
  d
  }




#' psiPrelec
#' 
#' This function is a transformation function that takes the log of the log of the inverse (which
#' can be used to generate the Prelec functions).
#' By design, it is robust to the inclusion of the occasional 0,
#' which it maps to a very small value, which can optionally be specified
#' It's appropriate for mapping (0, 1):-> (-inf, inf).
#' @param stimuli a vector of positive stimuli, between 0 and inf
#' @param smallValue A value to set nominal 0's to, to avoid errors in plotting simulated data, or aberrant responses
#' @return a vector containing "warped" stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiPrelec(1:100/100)
#' (0:1000/1000) %>% psiPrelec() %>% vanillaBayes() %>% psiPrelecInverse()  # Implements Prelec, 1998
psiPrelec <- function (stimuli, smallValue=10^-30) {
  UseMethod("psiPrelec", stimuli)
}

#' @export
psiPrelec.default <- function(stimuli, smallValue=10^-30){stimuli}

#' @export
psiPrelec.subjectiveLogLikelihood <- function(stimuli, smallValue){stimuli}


#' @export
psiPrelec.list <- function(stimuli, smallValue){
  mapply(psiPrelec, stimuli, smallValue, SIMPLIFY=FALSE)
}

#' @export
psiPrelec.numeric <- function(stimuli, smallValue=10^-30){
  stimuli[stimuli==1] <- 1-smallValue
  stimuli[stimuli==0] <- smallValue
  d <- log(log(1/stimuli))
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
#' @param smallValue A value to set nominal 0's to, to avoid errors in plotting simulated data, or aberrant responses
#' @return a vector containing "unwarped" stimuli identical to the warped stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentity
#' @export
#' @examples
#' psiIdentityInverse(-10:10)
psiIdentityInverse <- function (warpedStimuli, smallValue=10^-30) {
UseMethod("psiIdentityInverse", warpedStimuli)
}

#' @export
psiIdentityInverse.list <- function(warpedStimuli, smallValue=10^-30){
  mapply(psiIdentityInverse, warpedStimuli, SIMPLIFY=FALSE)
}

#' @export
psiIdentityInverse.default <- function(warpedStimuli, smallValue=10^-30){warpedStimuli}


#' @export
psiIdentityInverse.subjectiveLogLikelihood <- function(warpedStimuli, smallValue=10^-30){warpedStimuli}


#' @export
psiIdentityInverse.numeric <- function(warpedStimuli, smallValue=10^-30){
  warpedStimuli
}



#' psiLinearInverse
#' 
#' Creates mappings from -inf:+inf -> -inf:+inf, with a linear transformation mechanism
#' @param warpedStimuli a vector of stimuli, between -inf and inf
#' @param shift A scalar by which to (positively) offset the scaled values, or a list or vector of scalars
#' @param scaling A scalar multiplier to the (unshifted) values, or a list or vector of scalars
#' @return A vector containing unwarped stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiIdentity, psiLogOdds, psiPrelec, psiLog, psiLinear
#' @export
#' @examples
#' psiLinear(c(0.1, 0.2, 0.3))
#' psiLinear(-10:10)
#' psiLinear(-10:10, shift=10, scaling=2)
#' -10:10 %>% psiLinear(shift=2, scaling=0.5)
psiLinearInverse <- function (warpedStimuli, shift=0, scaling=1) {
  UseMethod("psiLinearInverse", warpedStimuli)
}

#' @export
psiLinearInverse.list <- function(warpedStimuli, shift=0, scaling=1){
  mapply(psiLinearInverse, warpedStimuli, shift, scaling, SIMPLIFY=FALSE)
}

#' @export
psiLogOdds.default <- function(stimuli, smallValue){stimuli}

#' @export
psiLinearInverse.subjectiveLogLikelihood <- function(warpedStimuli, shift=0, scaling=1){warpedStimuli}

#' @export
psiLinearInverse.numeric <- function(warpedStimuli, shift=0, scaling=1){
  (warpedStimuli-shift)/scaling 
}



#' psiLogInverse
#' 
#' This function is a transformation function that inverts the log...
#' It's literally just exp(stimuli).
#' It's appropriate for mapping  (-inf, inf):-> (0, +inf)
#' @param warpedStimuli a vector of stimuli, between -inf and inf
#' @param smallValue values smaller than this will be treated as this, to avoid anomalies near 0
#' @return A vector containing unwarped stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiLogInverse(psiLog(1:100)) # returns 1:100
#' psiLog(1:100) %>% psiLogInverse()
#' 1:1000 %>% psiLog() %>% vanillaBayes() %>% psiLogInverse()  # Implements Stevens Power Law
psiLogInverse <- function (warpedStimuli, smallValue=10^-30) {
  UseMethod("psiLogInverse", warpedStimuli )
}

#' @export
psiLogInverse.list <- function(warpedStimuli, smallValue=10^-30){
  mapply(psiLogInverse, warpedStimuli, smallValue, SIMPLIFY=FALSE)
}

#' @export
psiLogInverse.default <- function(warpedStimuli, smallValue=10^-30){warpedStimuli}


#' @export
psiLogInverse.subjectiveLogLikelihood <- function(warpedStimuli, smallValue=10^-30){warpedStimuli}

#' @export
psiLogInverse.numeric <- function(warpedStimuli, smallValue=10^-30){
  exp(warpedStimuli)
}

#' psiLogOddsInverse
#' 
#' This function is a transformation function that inverts the log odds.

#' It's appropriate for mapping (-inf, inf) :->  (0, 1)
#' As in all functions in this package, the parameters are set so that applying the same parameters to the main function 
#' and the inverse yields an identity.
#' @param warpedStimuli a vector of stimuli, between -inf and inf
#' @param smallValue A value to set nominal 0's to, to avoid errors in plotting simulated data, or aberrant responses. Accepted here to give the function and its inverse the same parameter set
#' @return A vector containing unwarped stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLog, psiIdentityInverse
#' @export
#' @examples
#' psiLogOdds(1:100/100)
#' (0:1000/1000) %>% psiLogOdds() %>% 
#'     vanillaBayes() %>% 
#'     psiLogOddsInverse()  # Implements Gonzales & Wu, 1996
psiLogOddsInverse <- function (warpedStimuli, smallValue=10^-30) {
  UseMethod("psiLogOddsInverse", warpedStimuli)
}

#' @export
psiLogOddsInverse.default <- function(warpedStimuli, smallValue=10^-30){warpedStimuli}

#' @export
psiLogOddsInverse.list <- function(warpedStimuli, smallValue=10^-30){
  mapply(psiLogOddsInverse, warpedStimuli, smallValue, SIMPLIFY=FALSE)
}

#' @export
psiLogOddsInverse.subjectiveLogLikelihood <- function(warpedStimuli, smallValue=10^-30){
  warpedStimuli
  }

#' @export
psiLogOddsInverse.numeric <- function(warpedStimuli, smallValue=10^-30){
  e <- exp(warpedStimuli) 
  d <- e/(1+e)
  d
}



#' psiPrelecInverse
#' 
#' This function is a transformation function inverts the Prelec transformation, ln(ln(1/s))
#' By design, it is robust to the inclusion of the occasional 0,
#' which it maps to a very small value, which can optionally be specified
#' It's appropriate for mapping (0, 1):-> (-inf, inf).
#' @param warpedStimuli a vector of stimuli, between -inf and inf
#' @param smallValue A value to set nominal 0's to, to avoid errors in plotting simulated data, or aberrant responses. Accepted here to give the function and its inverse the same parameter set
#' @return A vector containing unwarped stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiLinear, psiLogOdds, psiPrelec, psiLogInverse, psiIdentityInverse
#' @export
#' @examples
#' psiPrelecInverse(-10:10)
#' (0:1000/1000) %>% psiPrelec() %>% 
#'     vanillaBayes() %>% 
#'     psiPrelecInverse()  # Implements Prelec, 1998
psiPrelecInverse <- function (warpedStimuli, smallValue=10^-30) {
  UseMethod("psiPrelecInverse", warpedStimuli)
}

#' @export
psiPrelecInverse.default <- function(warpedStimuli, smallValue=10^-30){warpedStimuli}

#' @export
psiPrelecInverse.list <- function(warpedStimuli, smallValue){
  mapply(psiPrelec, warpedStimuli, smallValue, SIMPLIFY=FALSE)
}

#' @export
psiPrelecInverse.subjectiveLogLikelihood <- function(warpedStimuli, smallValue=10^-30){warpedStimuli}

#' @export
psiPrelecInverse.numeric <- function(warpedStimuli, smallValue=10^-30){
  exp(-exp(warpedStimuli))
  
}


# MULTICYCLING

# Multicycling is an advanced art.  Basically, we 'multicycle' by a pair of functions and inverses 
# that carve spaces into regions that are already bounded at lynchpin points.  These can include 'inf' and '-inf'


#' multiCycle
#' 
#' Divides a mapping, of any sort (bounded, semi-bounded, or unbounded) into a finite number of distinct regions,
#' by creating 'references' inside the space. Functionally, multiCycle simply divides a vector of stimuli into 
#' A list of vectors of values based on their proportional distances between two references.  Usually, these
#' are then fed into a psi function which creates a multiCycle space containing multiple infinite unbounded regions
#' which are conceptually adjacent
#'
#' For instance, a spatial region |----------------|
#' might be divided by a line of vertical symmetry into two regions
#'  |--------||--------|
#'  And each be turned int an unbounded Prelec region
#'  <--------><-------->
#'  with a command like -10:10 %>% multicycle(0) %>% psiPrelec
#' @param stimuli a vector of stimuli, between -inf and inf, or a list of vectors of stimuli
#' @param references A vector of values.  This should include -Inf and Inf, but these are implicitly added
#' if they are omitted
#' @return A vector containing warped stimuli, or a list of vectors
#' @keywords psi perceptual tranformations
#' @seealso psiIdentity, multiCycleInverse
#' @export
#' @examples
#' multiCycle(-99:100, c(-100, 0, 100))
#' (-99:100/100) %>% multiCycle(references= c(-1, 0, 1)) %>% 
#'     psiLogOdds() %>% vanillaBayes() %>% 
#'     psiLogOddsInverse() %>% multiCycleInverse(references=c(-1, 0, 1)) 
#'     # Implements Landy et al's model of one-dimensional spatial memory, with fixed boundaries
#' plot(-99:100, unlist(multiCycle(-99:100, c(-100, -50, 0, 50, 100))))
#' plot(-99:100/100, (-99:100/100) %>% multiCycle(references= c(-1, 0, 1)) %>% 
#'     psiLogOdds() %>% vanillaBayes(kappa=c(-.8, .8)) %>% 
#'     psiLogOddsInverse() %>%
#'      multiCycleInverse(references=c(-1, 0, 1))-(-99:100/100),
#'           ylab="bias", xlab="stimulus");abline(0,0)
multiCycle <- function (stimuli, references=c(0)) {
  UseMethod("multiCycle", stimuli)
}

#' @export
multiCycle.list <- function(stimuli, references=c(0)){
  mapply(multiCycle, stimuli, references, SIMPLIFY=FALSE)
}

#' @export
multiCycle.subjectiveLogLikelihood <- function(stimuli, references=c(0)){stimuli}

#' @export
multiCycle.default <- function(stimuli, references=c(0)){stimuli}


#' @export
multiCycle.numeric <- function(stimuli, references=c(0)){
  multiCycleScalingFunction <- function(stimuli, left, right){
    if(left==-Inf && right==Inf){
      stimuli
    } else if(left == -Inf){
      stimuli - right
    } else if(right == Inf){
      stimuli - left
    } else {
      (stimuli-left)/(right-left)    
    }
  }
  
  if(references[1]!= -Inf){references <- c(-Inf, references)}
  if(references[length(references)]!= Inf){references <- c(references, Inf)}
  brokenVersion <- split(stimuli, cut(stimuli,breaks=references, ordered_result=TRUE))
  mapply(multiCycleScalingFunction,
         brokenVersion
         , references[1:length(references)-1]
         , references[2:length(references)]
         , SIMPLIFY=FALSE)
}  



#' multiCycleInverse
#' 
#' recombines a mapping, of any sort (bounded, semi-bounded, or unbounded) into single region,
#' by eliminating 'references' inside the space. Functionally, multiCycleInvers simply undoes a vector of stimuli into 
#'
#' For instance, a spatial region |----------------|
#' might be divided by a line of vertical symmetry into two regions
#'  |--------||--------|
#'  And each be turned int an unbounded Prelec region
#'  <--------><-------->
#'  with a command like -10:10 %>% multicycle(0) %>% psiPrelec
#'  multiCycleInverse would then recombine these values into a single range.
#'
#'  KEY LIMITATION: Right now, multiCycleInverse can only 'handle' one layer of recursion. That is, you can't make calls like
#'  -10:10 %>% multiCycle %>% multiCycle %> multiCycleInverse %>% multiCycleInverse
#'  . Sorry.
#' @param warpedStimuli a vector of stimuli, between -inf and inf, or a list of vectors of stimuli
#' @param references A vector of values.  This should include -Inf and Inf, but these are implicitly added
#' if they are omitted
#' @return A vector containing warped stimuli, or a list of vectors
#' @keywords psi perceptual tranformations
#' @seealso psiIdentity, psiLogOdds, psiPrelec, psiLog, psiLinearInverse
#' @export
#' @examples
#' multiCycle(-99:100, c(-100, 0, 100))
#' plot(-99:100, unlist(multiCycle(-99:100, c(-100, -50, 0, 50, 100))))
#' (-99:100/100) %>% multiCycle(-1, 0, 1) %>% 
#'     psiLogOdds() %>% vanillaBayes() %>% 
#'     psiLogOddsInverse() %>% 
#'     multiCycleInverse(-1, 0, 1) # Implements Landy et al's model 
#'     # of one-dimensional spatial memory, with fixed boundaries
multiCycleInverse <- function (warpedStimuli, references=c(0)) {
  UseMethod("multiCycleInverse", warpedStimuli)
}

#' @export
multiCycleInverse.numeric <- function(warpedStimuli, references=c(0)){sum(warpedStimuli)}



#' @export
multiCycleInverse.subjectiveLogLikelihood <- function(warpedStimuli, references=c(0)){
  warpedStimuli
  }

#' @export
multiCycleInverse.list <- function(warpedStimuli, references=c(0)){
  multiCycleInverseScalingFunction <- function(warpedStimuli, left, right){
    if("subjectiveLogLikelihood" %in% class(warpedStimuli)){return(sum(warpedStimuli))}
    if(left==-Inf && right==Inf){
      warpedStimuli
    } else if(left == -Inf){
      warpedStimuli + right
    } else if(right == Inf){
      warpedStimuli + left
    } else {
      warpedStimuli*(right-left) + left
    }
  }
  if(references[1]!= -Inf){references <- c(-Inf, references)}
  if(references[length(references)]!= Inf){references <- c(references, Inf)}
  unlist(mapply(multiCycleInverseScalingFunction,
         warpedStimuli 
         , references[1:length(references)-1]
         , references[2:length(references)]
         , SIMPLIFY=TRUE)
         , recursive=FALSE
         , use.names=FALSE)
}  






# 

#' unicycle
#' 
#' Divides a mapping, of any sort (bounded, semi-bounded, or unbounded) into a finite number of distinct regions,
#' by creating 'references' inside the space. Functionally, unicycle simply divides a vector of stimuli into 
#' A list of vectors of values based on their proportional distances between two references.  Usually, these
#' are then fed into a psi function which creates a unbounded space containing multiple infinite unbounded regions
#' which are conceptually adjacent.
#' 
#' Unicycling is a simplification of multicycling. Unicycling (and its inverse), like multicycling,
#  carves spaces into a region that is bounded at lynchpin points.  These can include 'inf' and '-inf'.  
# The difference between unicycle and multicycle is that while multicycle divides regions up into 
# similar kinds of adjacent spaces, unicycle puts 0 probability density in regions outside the focal region. 

#'
#' For instance, a spatial region |----------------|
#' might be divided into three regions
#'  |=======||--------||=========|
#'  Where the central |--------| region is turned into a preLec function, and the 
#'  regions outside are 'no op' regions, in whihc no predictions can be made, and no probability density is placed 
#'  (or, to keep things from underflowying, very very little).
#'  with a command like -10:10 %>% uniCycle(0) %>% psiPrelec
#' @param stimuli a vector of stimuli, between -inf and inf, or a list of vectors of stimuli
#' @param references A vector of values.  This should include -Inf and Inf, but these are implicitly added
#' if they are omitted
#' @return A vector containing warped stimuli, or a list of vectors
#' @keywords psi perceptual tranformations
#' @seealso psiIdentity, multiCycle
#' @export
#' @examples
#' uniCycle(-99:100, c(-100, 0, 100))
#' (-99:100/100) %>% multiCycle(references= c(-1, 0, 1)) %>% 
#'     psiLogOdds() %>% vanillaBayes() %>% 
#'     psiLogOddsInverse() %>% multiCycleInverse(references=c(-1, 0, 1)) 
#'     # Implements Landy et al's model of one-dimensional spatial memory, with fixed boundaries
#' plot(-99:100, unlist(multiCycle(-99:100, c(-100, -50, 0, 50, 100))))
#' plot(-99:100/100, (-99:100/100) %>% multiCycle(references= c(-1, 0, 1)) %>% 
#'     psiLogOdds() %>% vanillaBayes(kappa=c(-.8, .8)) %>% 
#'     psiLogOddsInverse() %>%
#'      uniCycleInverse(references=c(-1, 0, 1))-(-99:100/100),
#'           ylab="bias", xlab="stimulus");abline(0,0)
uniCycle <- function (stimuli, references=c(0)) {
  UseMethod("uniCycle", stimuli)
}

#' @export
uniCycle.list <- function(stimuli, references=c(0)){
  mapply(uniCycle, stimuli, references, SIMPLIFY=FALSE)
}

#' @export
uniCycle.subjectiveLogLikelihood <- function(stimuli, references=c(0)){stimuli}

#' @export
uniCycle.default <- function(stimuli, references=c(0)){stimuli}


#' @export
uniCycle.numeric <- function(stimuli, references=c(0)){
  uniCycleScalingFunction <- function(stimuli, left, right){
    if(left==-Inf && right==Inf){
      stimuli
    } else if(left == -Inf){
      stimuli - right
    } else if(right == Inf){
      stimuli - left
    } else {
      (stimuli-left)/(right-left)    
    }
  }
  
  if(references[1]!= -Inf){references <- c(-Inf, references)}
  if(references[length(references)]!= Inf){references <- c(references, Inf)}
  brokenVersion <- split(stimuli, cut(stimuli,breaks=references, ordered_result=TRUE))
  mapply(uniCycleScalingFunction,
         brokenVersion
         , references[1:length(references)-1]
         , references[2:length(references)]
         , SIMPLIFY=FALSE)
}  



#' uniCycleInverse
#' 
#' recombines a mapping, of any sort (bounded, semi-bounded, or unbounded) into single region,
#' by eliminating 'references' inside the space. Functionally, uniCycleInverse simply undoes a vector of stimuli into 
#'
#' For instance, a spatial region |----------------|
#' might be divided into three regions
#'  |=======||--------||=========|

#'  And the middle be turned into an unbounded Prelec region
#'  |=======|<-------->|=========|
#'  with a command like -10:10 %>% unicycle(0) %>% psiPrelec
#'  uniCycleInverse would then recombine these values into a single range.
#'
#'  KEY LIMITATION: Right now, uniCycleInverse can only 'handle' one layer of recursion. That is, you can't make calls like
#'  -10:10 %>% uniCycle %>% uniCycle %> uniCycleInverse %>% uniCycleInverse
#'  . Sorry.
#' @param warpedStimuli a vector of stimuli, between -inf and inf, or a list of vectors of stimuli
#' @param references A vector of values.  This should include -Inf and Inf, but these are implicitly added
#' if they are omitted
#' @return A vector containing warped stimuli, or a list of vectors
#' @keywords psi perceptual tranformations
#' @seealso psiIdentity, psiLogOdds, psiPrelec, psiLog, psiLinearInverse
#' @export
#' @examples
#' uniCycle(-99:100, c(-100, 0, 100))
#' plot(-99:100, unlist(uniCycle(-99:100, c(-100, -50, 0, 50, 100))))
#' (-99:100/100) %>% uniCycle(-1, 0, 1) %>% 
#'     psiLogOdds() %>% vanillaBayes() %>% 
#'     psiLogOddsInverse() %>% 
#'     uniCycleInverse(-1, 0, 1) # Implements Landy et al's model 
#'     # of one-dimensional spatial memory, with fixed boundaries
uniCycleInverse <- function (warpedStimuli, references=c(0)) {
  UseMethod("uniCycleInverse", warpedStimuli)
}

#' @export
uniCycleInverse.numeric <- function(warpedStimuli, references=c(0)){sum(warpedStimuli)}



#' @export
uniCycleInverse.subjectiveLogLikelihood <- function(warpedStimuli, references=c(0)){
  warpedStimuli
}

#' @export
uniCycleInverse.list <- function(warpedStimuli, references=c(0)){
  uniCycleInverseScalingFunction <- function(warpedStimuli, left, right){
    if("subjectiveLogLikelihood" %in% class(warpedStimuli)){return(sum(warpedStimuli))}
    if(left==-Inf && right==Inf){
      warpedStimuli
    } else if(left == -Inf){
      warpedStimuli + right
    } else if(right == Inf){
      warpedStimuli + left
    } else {
      warpedStimuli*(right-left) + left
    }
  }
  if(references[1]!= -Inf){references <- c(-Inf, references)}
  if(references[length(references)]!= Inf){references <- c(references, Inf)}
  unlist(mapply(uniCycleInverseScalingFunction,
                warpedStimuli 
                , references[1:length(references)-1]
                , references[2:length(references)]
                , SIMPLIFY=TRUE)
         , recursive=FALSE
         , use.names=FALSE)
}  



