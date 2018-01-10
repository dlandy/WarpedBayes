# classicPsychophysics models



#' classicStevensPowerLaw
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param exponent the power to raise the stimulus to, equivalent to the ratio of taus in the reconsturcted version 
#' @param proportionality the constant multiple on the stimulus value
#' @param responses an optional vector of responses.
#' @param mode a mode. can validly be "prediction", or "logLikelihoodOfResponses". If the latter, responses must be given, and tau must be given
#' @param tau if log-likelihood is to be calculated, we must decide the predicted variability of responses. This is captured in tau, the precision
#' @return A vector the transformed stimuli 
#' @seealso classicGonzalezWu, classicSpence
#' @export
#' @examples
#' classicStevensPowerLaw(1:100)
#' classicStevensPowerLaw(1:100, exponent=0.5, proportionality=2)
#' classicStevensPowerLaw(1:100, exponent=0.9, proportionality=2, responses=2*(1:100)^.9)
classicStevensPowerLaw <- function(stimuli
                                   , exponent = 1
                                   , proportionality = 1
                                   , responses = NULL
                                   , tau = 1
                                   , mode="prediction"){
  predictions <- proportionality * stimuli^exponent
  if(mode=="prediction"){
    return(predictions)
  } else{ 
    return(0-sum(log(dnorm(responses-predictions, sd=1/sqrt(tau)))))
  }
}


