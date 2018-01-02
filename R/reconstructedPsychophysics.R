# Reconstructed classical functions, using Bayesian parameterizations


#' bayesianStevensPowerLaw
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param kappa The location of the category
#' @param tauStimuli The precision of the stimulus traces: may be a single number or a vector
#' @param tauCategory The precision of the category distribution
#' @param responses an optional vector of responses. If responses are given, the return value is the 

#' @return A vector the transformed stimuli
#' @seealso bayesianHuttenlocherSpatialMemory
#' @export
#' @examples
#' bayesianStevensPowerLaw(1:100)
#' bayesianStevensPowerLaw(1:100, kappa=1, tauStimuli=2)
#' bayesianStevensPowerLaw(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianStevensPowerLaw <- function(stimuli,  kappa=0, tauStimuli=1, tauCategory=1, responses="none"){
    stimuli %>% psiLog %>% vanillaBayes(kappa=kappa
                                        , tauStimuli=tauStimuli
                                        , tauCategory=tauCategory
                                        , responses=psiLog(responses)) %>% psiLogInverse()
}


#' bayesianSpatialMemoryHuttenlocher
#' 
#' A simple package that assembles the symmetric model used by Huttenlocher and colleagues analyze spatial 
#' estimations from memory in one dimension
#' @param stimuli a vector of stimuli, between -inf and inf
#' @param kappa The location of the categories (presumed symmetric on both sides around 0)
#' @param tauStimuli The precision of the stimulus traces: should be a single number
#' @param tauCategory The precision of the category distribution: should be a single number
#' @param boundaries The subject-specific location of the boundaries: may bear any relation to true stimuli, except that it 
#' should not leave real data outside the boundaries
#' @param responses an optional vector of responses. If responses are given, the return value is the logLikelihood of the responses given the parameters
#' @return A vector the transformed stimuli, or the logLikelihood of them.
#' @seealso psiIdentity, multiCycleInverse
#' @export
#' @examples
#' bayesianSpatialMemoryHuttenlocher(-99:100/100)
#' bayesianSpatialMemoryHuttenlocher(-99:100/100, kappa=1, tauStimuli=2)
#' bayesianSpatialMemoryHuttenlocher(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianSpatialMemoryHuttenlocher <- function(stimuli
                                              , kappa=0.5
                                              , tauStimuli=1
                                              , tauCategory=1
                                              , boundary = 1
                                              , leftBoundary = -1*boundary
                                              , rightBoundary = boundary
                                              , center = 0
                                              , responses="none"){
  refs <- c(leftBoundary, center, rightBoundary)
  stimuli %>% multiCycle(references = refs) %>%
    vanillaBayes(kappa=kappa
                 , tauStimuli=tauStimuli
                 , tauCategory=tauCategory
                 , responses=multiCycle(responses, references = refs))  %>% 
    multiCycleInverse(references = refs)
}




#' bayesianSpatialMemoryLandyCrawfordCorbin2017
#' 
#' A simple package that assembles the symmetric model used by Huttenlocher and colleagues analyze spatial 
#' estimations from memory in one dimension
#' @param stimuli a vector of stimuli, between -inf and inf
#' @param kappa The location of the categories (presumed symmetric on both sides around 0)
#' @param tauStimuli The precision of the stimulus traces: should be a single number
#' @param tauCategory The precision of the category distribution: should be a single number
#' @param boundaries The subject-specific location of the boundaries: may bear any relation to true stimuli, except that it 
#' should not leave real data outside the boundaries
#' @param responses an optional vector of responses. If responses are given, the return value is the logLikelihood of the responses given the parameters
#' @return A vector the transformed stimuli, or the logLikelihood of them.
#' @seealso psiIdentity, multiCycleInverse
#' @export
#' @examples
#' bayesianSpatialMemoryLandyCrawfordCorbin2017(-99:100/100)
#' bayesianSpatialMemoryLandyCrawfordCorbin2017(-99:100/100, kappa=1, tauStimuli=2)
#' bayesianSpatialMemoryLandyCrawfordCorbin2017(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianSpatialMemoryLandyCrawfordCorbin2017 <- function(stimuli
                                              , kappa=0.5
                                              , tauStimuli=1
                                              , tauCategory=1
                                              , boundary = 1
                                              , leftBoundary = -1*boundary
                                              , rightBoundary = boundary
                                              , center = 0
                                              , responses="none"){
  refs <- c(leftBoundary, center, rightBoundary)
  kapp <- c(kappa, 1-kappa)
  stimuli %>% multiCycle(references = refs) %>%
    psiLogOdds() %>% 
    vanillaBayes(kappa=kappa
                 , tauStimuli=tauStimuli
                 , tauCategory=tauCategory
                 , responses=multiCycle(responses, references = refs))  %>% 
    psiLogOddsInverse() %>% 
    multiCycleInverse(references = refs)
}




#Helper functions that fit data!
#' @export
fitWarpedBayesModel <- function(model, stimuli, responses
                , initialPars
                , parNames=initialPars
                , control=list(maxit=5000, reltol = 10e-120)
){
  require(tidyverse)
  fitFunction <- function(pars){
    do.call(model, append(append(list(stimuli=stimuli), pars), list(responses=responses)))
  }
  result <- optim(initialPars, fitFunction, control=control )
  predictions <- do.call(model, append(append(list(stimuli=stimuli), result$par), list(responses="simulation")))
  a <- tibble(
    stimulus = stimuli
    , response = responses
    , prediction = predictions
    , value = result$value
    , counts=result$counts[1]
    , convergence = result$convergence
  )
  for(i in 1:length(result$par)){
    a[names(result$par[i])] <- result$par[i]
  }
  
  
  a
}



