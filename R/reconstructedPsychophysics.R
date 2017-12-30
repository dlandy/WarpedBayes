# Reconstructed classical functions, using Bayesian parameterizations


#' bayesianStevensPowerLaw
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param kappa The location of the category
#' @param tauStimuli The precision of the stimulus traces: may be a single number or a vector
#' @param tauCategory The precision of the category distribution
#' @param responses an optional vector of responses. If responses are given, the return value is the 

#' @return A vector the transformed stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiIdentity, multiCycleInverse
#' @export
#' @examples
#' bayesianStevensPowerLaw(1:100)
#' bayesianStevensPowerLaw(1:100, kappa=1, tauStimuli=2)
#' bayesianStevensPowerLaw(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianStevensPowerLaw <- function(stimuli,  kappa=0, tauStimuli=1, tauCategory=1, responses="none"){
  result <- stimuli %>% psiLog %>% vanillaBayes(kappa=kappa, tauStimuli=tauStimuli, tauCategory=tauCategory, responses=responses) 
  if(length(responses)==1 && responses=="none"){
   result %>% psiLogInverse
  } else {
    result
  }
}