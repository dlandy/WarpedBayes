# classicPsychophysics models



#' stevensPowerLaw
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param beta the power to raise the stimulus to, equivalent to the ratio of taus in the reconstucted version 
#' @param delta proportionality constant--the constant multiple on the stimulus value
#' @param responses an optional vector of responses.
#' @param mode a mode. can validly be "prediction", or "logLikelihoodOfResponses". If the latter, responses must be given, and tau must be given
#' @param tau if log-likelihood is to be calculated, we must decide the predicted variability of responses. This is captured in tau, the precision
#' @return A vector the transformed stimuli 
#' @seealso gonzalezWu, spence
#' @export
#' @examples
#' stevensPowerLaw(1:100)
#' stevensPowerLaw(1:100, beta=0.5, delta=2)
#' stevensPowerLaw(1:100, beta=0.9, delta=2, responses=2*(1:100)^.9, mode="logLikelihoodOfResponses)
stevensPowerLaw <- function(stimuli
                                   , beta = 1
                                   , delta = 1
                                   , responses = NULL
                                   , tau = 1
                                   , mode="prediction"){
  predictions <- delta * stimuli^beta
  if(mode=="prediction"){
    return(predictions)
  } else{ 
    return(0-sum(log(dnorm(responses-predictions, sd=1/sqrt(tau)))))
  }
}







#' spence
#' @param stimuli a vector of stimuli, between 0 and 1
#' @param beta the power to raise the stimulus to, equivalent to the ratio of taus in the reconsturcted version 
#' @param responses an optional vector of responses.
#' @param mode a mode. can validly be "prediction", or "logLikelihoodOfResponses". If the latter, responses must be given, and tau must be given
#' @param tau if log-likelihood is to be calculated, we must decide the predicted variability of responses. This is captured in tau, the precision
#' @return A vector the transformed stimuli 
#' @seealso classicGonzalezWu, spence
#' @export
#' @examples
#' spence(seq(0.00, 1, 0.01))
#' spence(seq(0.00, 1, 0.01), beta=0.5)
#' spence(seq(0.00, 1, 0.01), beta=0.9, responses=spence(seq(0.00, 1, 0.01), 0.8), mode="logLikelihoodOfResponses)
spence <- function(stimuli
                                   , beta = 1
                                   , responses = NULL
                                   , tau = 1
                                   , mode="prediction"){
  predictions <- stimuli^beta/(stimuli^beta + (1-stimuli)^beta)
  if(mode=="prediction"){
    return(predictions)
  } else{ 
    return(0-sum(log(dnorm(responses-predictions, sd=1/sqrt(tau)))))
  }
}


#' gonzalezWu
#' @param stimuli a vector of stimuli, between 0 and 1
#' @param beta the power to raise the stimulus to, equivalent to the ratio of taus in the reconsturcted version 
#' @param delta the power to raise the stimulus to, equivalent to the ratio of taus in the reconsturcted version 
#' @param responses an optional vector of responses.
#' @param mode a mode. can validly be "prediction", or "logLikelihoodOfResponses". If the latter, responses must be given, and tau must be given
#' @param tau if log-likelihood is to be calculated, we must decide the predicted variability of responses. This is captured in tau, the precision
#' @return A vector the transformed stimuli 
#' @seealso classicGonzalezWu, spence
#' @export
#' @examples
#' gonzalezWu(seq(0.00, 1, 0.01))
#' gonzalezWu(seq(0.00, 1, 0.01), beta=0.5)
#' gonzalezWu(seq(0.00, 1, 0.01), beta=0.9, responses=gonzalezWu(seq(0.00, 1, 0.01), 0.8), mode="logLikelihoodOfResponses)
gonzalezWu <- function(stimuli
                   , beta = 1
                   , delta = 1
                   , responses = NULL
                   , tau = 1
                   , mode="prediction"){
  predictions <- delta*stimuli^beta/(delta*stimuli^beta + (1-stimuli)^beta)
  if(mode=="prediction"){
    return(predictions)
  } else{ 
    return(0-sum(log(dnorm(responses-predictions, sd=1/sqrt(tau)))))
  }
}


#' landyBrowerBiasedSpence
#' @param stimuli a vector of stimuli, between 0 and 1
#' @param beta the power to raise the stimulus to, equivalent to the ratio of taus in the reconsturcted version 
#' @param bias a fixed shift (in log odds space) of the stimuli
#' @param responses an optional vector of responses.
#' @param mode a mode. can validly be "prediction", or "logLikelihoodOfResponses". If the latter, responses must be given, and tau must be given
#' @param tau if log-likelihood is to be calculated, we must decide the predicted variability of responses. This is captured in tau, the precision
#' @return A vector the transformed stimuli 
#' @seealso classicGonzalezWu, spence
#' @export
#' @examples
#' landyBrowerBiasedSpence(seq(0.00, 1, 0.01))
#' landyBrowerBiasedSpence(seq(0.00, 1, 0.01), beta=0.5)
#' landyBrowerBiasedSpence(seq(0.00, 1, 0.01), beta=0.9, bias=1, responses=gonzalezWu(seq(0.00, 1, 0.01), 0.8), mode="logLikelihoodOfResponses)
landyBrowerBiasedSpence <- function(stimuli
                       , beta = 1
                       , bias = 0
                       , responses = NULL
                       , tau = 1
                       , mode="prediction"){
  predictions <- bias^beta*stimuli^beta/(bias^beta*stimuli^beta + (1-stimuli)^beta)
  if(mode=="prediction"){
    return(predictions)
  } else{ 
    return(0-sum(log(dnorm(responses-predictions, sd=1/sqrt(tau)))))
  }
}





