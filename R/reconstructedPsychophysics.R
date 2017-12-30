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


#' bayesianHuttenlocherSpatialMemory
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
#' bayesianHuttenlocherSpatialMemory(-99:100/100)
#' bayesianHuttenlocherSpatialMemory(-99:100/100, kappa=1, tauStimuli=2)
#' bayesianHuttenlocherSpatialMemory(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianHuttenlocherSpatialMemory <- function(stimuli
                                              ,  kappa=0
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





#' bayesianSpatialEstimationLandyCrawfordCorbin2017
#' 
#' A simple package that assembles the symmetric model used by Landy et al 2017 to analyze spatial 
#' estimations from memory in one dimension
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param kappa The location of the category
#' @param tauStimuli The precision of the stimulus traces: may be a single number or a vector
#' @param tauCategory The precision of the category distribution
#' @param responses an optional vector of responses. If responses are given, the return value is the 

#' @return A vector the transformed stimuli
#' @keywords psi perceptual tranformations
#' @seealso bayesianStevensPowerLaw
#' @export
#' @examples
#' bayesianStevensPowerLaw(1:100)
#' bayesianStevensPowerLaw(1:100, kappa=1, tauStimuli=2)
#' bayesianStevensPowerLaw(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianSpatialEstimationLandyCrawfordCorbin2017 <- function(stimuli,  kappa=0, tauStimuli=1, tauCategory=1, responses="none"){
  stimuli %>% psiLog %>% vanillaBayes(kappa=kappa
                                      , tauStimuli=tauStimuli
                                      , tauCategory=tauCategory
                                      , responses=psiLog(responses)) %>% psiLogInverse()
}




