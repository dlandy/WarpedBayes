# Reconstructed classical functions, using Bayesian parameterizations


#' bayesianStevensPowerLaw
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param kappa The location of the category
#' @param kappaObjective An alternative specification giving kappa in objective units
#' @param tauStimuli The precision of the stimulus traces: may be a single number or a vector
#' @param tauCategory The precision of the category distribution
#' @param responses an optional vector of responses.
#' If responses are given, the return value is the logLikelihood of the responses given the parameters

#' @return A vector the transformed stimuli
#' @seealso bayesianHuttenlocherSpatialMemory
#' @export
#' @examples
#' bayesianStevensPowerLaw(1:100)
#' bayesianStevensPowerLaw(1:100, kappa=1, tauStimuli=2)
#' bayesianStevensPowerLaw(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianStevensPowerLaw <- function(stimuli
                                    , kappaObjective = 1
                                    , kappa=psiLog(kappaObjective), tauStimuli=1, tauCategory=1, responses="prediction"){
    stimuli %>% psiLog %>% vanillaBayes(kappa=kappa
                                        , tauStimuli=tauStimuli
                                        , tauCategory=tauCategory
                                        , responses=psiLog(responses)) %>% psiLogInverse()
}


#' bayesianSpatialMemoryHuttenlocher
#' 
#' A simple package that assembles the symmetric model used by Huttenlocher and colleagues analyze spatial 
#' estimations from memory in one dimension.
#' @param stimuli a vector of stimuli, in public units between -inf and inf
#' @param kappa The location of the right-side category (presumed symmetric on both sides around the midline of the screen)
#' @param kappaObjective An alternative specification giving kappa in objective units
#' @param tauStimuli The precision of the stimulus traces: should be a single number
#' @param tauCategory The precision of the category distribution: should be a single number
#' @param boundary The subject-specific location of the boundaries: may bear any relation to true stimuli, except that it 
#' should not leave real data outside the boundaries
#' @param leftBoundary The location of the posited (or fitted) psychological left-hand boundary of the screen. Defaults to -1 * 'boundary'
#' @param rightBoundary The location of the posited (or fitted) psychological right-hand boundary of the screen. Defaults to 'boundary'
#' @param center The posited (or fitted) psychological center of the screen (in public units: should be near the true center)
#' @param responses an optional vector of responses. If responses are given, the return value is the logLikelihood of the responses given the parameters
#' @return A vector the transformed stimuli, or the logLikelihood of them.
#' @details This package 
#' @seealso psiIdentity, multiCycleInverse
#' @export
#' @examples
#' bayesianSpatialMemoryHuttenlocher(-99:100/100)
#' bayesianSpatialMemoryHuttenlocher(-99:100/100, kappa=1, tauStimuli=2)
#' bayesianSpatialMemoryHuttenlocher(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianSpatialMemoryHuttenlocher <- function(stimuli
                                              , kappaObjective = 0.5
                                              , kappa=psiLogOdds(kappaObjective)
                                              , tauStimuli=1
                                              , tauCategory=1
                                              , boundary = 1
                                              , leftBoundary = -1*boundary
                                              , rightBoundary = boundary
                                              , center = 0
                                              , responses="prediction"){
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
#' @param kappa The location of the categories (presumed symmetric on both sides around the midline of the screen)
#' @param kappaObjective The location of the categories measured in objective units
#' @param tauStimuli The precision of the stimulus traces: should be a single number
#' @param tauCategory The precision of the category distribution: should be a single number
#' @param boundary The subject-specific location of the boundaries: may bear any relation to true stimuli, except that it 
#' should not leave real data outside the boundaries
#' @param leftBoundary The location of the posited (or fitted) psychological left-hand boundary of the screen. Defaults to -1 * 'boundary'
#' @param rightBoundary The location of the posited (or fitted) psychological right-hand boundary of the screen. Defaults to 'boundary'
#' @param center The posited (or fitted) psychological center of the screen (in public units: should be near the true center)
#' @param responses an optional vector of responses. If responses are given, the return value is the logLikelihood of the responses given the parameters
#' @return A vector the transformed stimuli, or the logLikelihood of them.
#' @seealso psiIdentity, multiCycleInverse
#' @export
#' @examples
#' bayesianSpatialMemoryLandyCrawfordCorbin2017(-99:100/100)
#' bayesianSpatialMemoryLandyCrawfordCorbin2017(-99:100/100, kappa=1, tauStimuli=2)
#' bayesianSpatialMemoryLandyCrawfordCorbin2017(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianSpatialMemoryLandyCrawfordCorbin2017 <- function(stimuli
                                              , kappaObjective = 0.5
                                              , kappa=psiLogOdds(kappaObjective)
                                              , tauStimuli=1
                                              , tauCategory=1
                                              , boundary = 1
                                              , leftBoundary = -1*boundary
                                              , rightBoundary = boundary
                                              , center = 0
                                              , responses="prediction"){
  refs <- c(leftBoundary, center, rightBoundary)
  kappas <- c(kappa, 1-kappa)
  stimuli %>% multiCycle(references = refs) %>%
    psiLogOdds() %>% 
    vanillaBayes(kappa=kappas
                 , tauStimuli=tauStimuli
                 , tauCategory=tauCategory
                 , responses=multiCycle(responses, references = refs) %>% psiLogOdds()
                 )  %>% 
    psiLogOddsInverse() %>% 
    multiCycleInverse(references = refs)
}





#' bayesianGonzalezWu
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param kappa The location of the category in subjective space (from -inf to inf). Defaults to 0
#' @param kappaObjective An alternative specification for the kappa location, situated in objective measures.
#' @param tauStimuli The precision of the stimulus traces: may be a single number or a vector
#' @param tauCategory The precision of the category distribution
#' @param leftBoundaryObjective The location (in objective units) of the posited (or fitted) psychological left-hand boundary of legal 
#' stimulus values. Defaults to a small amount less than 0
#' @param rightBoundaryObjective The location (in objective units) of the posited (or fitted) psychological right-hand boundary of legal 
#' stimulus values. Defaults to a small amount more than 1
#' @param leftBoundaryExpansion How much beyond the minimium stimulus/response value should the left boundary be expanded?
#' @param rightBoundaryExpansion How much beyond the maximum stimulus/response value should the right boundary be expanded?
#' @param smallValue a small amount, by which to default boundary expansion
#' @param minValue The lowest range value (in objective units).  Should be at least as small as the smallest stimulus and response.
#' @param maxValue The highest range value  (in objective units).  Should be at least as large as the largest stimulus and response.
#' @param responses an optional vector of responses
#' @param mode What aspect should the function calculate? Legel choices include "prediction", "simulation", and "logLikelihoodOfResponses"
#' If responses are given, the return value is the logLikelihood of the responses given the parameters
#' @return A vector the transformed stimuli
#' @seealso bayesianHuttenlocherSpatialMemory
#' @export
#' @examples
#' bayesianGonzalezWu(1:100/100)
#' bayesianGonzalezWu(1:100/100, kappa=1, tauStimuli=2)
#' bayesianGonzalezWu(1:100/100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianGonzalezWu <- function(stimuli
                               , kappa = psiLogOdds(multiCycle(kappaObjective, references=c(leftBoundaryObjective, rightBoundaryObjective)))
                               , kappaObjective = 0.5
                               , tauStimuli = 1
                               , tauCategory = 1
                               , leftBoundaryObjective = 0-smallValue
                               , rightBoundaryObjective = 1 + smallValue
                               , rightBoundaryExpansion = log((rightBoundaryObjective-maxValue)/(maxValue - minValue))
                               , leftBoundaryExpansion  = log((minValue-leftBoundaryObjective )/(maxValue - minValue))
                               , minValue = min(c(stimuli, responses), na.rm=T)
                               , maxValue = max(c(stimuli, responses), na.rm=T)
                               , smallValue = 10^-10
                               , responses=NULL
                               , mode = "prediction"){
  # Safety check parameters
  if(minValue <= leftBoundaryObjective){
    warning("leftBoundaryObjective (", leftBoundaryObjective, ") larger than smallest stimulus (", minValue , ")")
    return(10^5+abs(minValue-leftBoundaryObjective)) # Return a large value for convenience for optim
  } else if(maxValue > rightBoundaryObjective){
    warning("rightBoundaryObjective (", rightBoundaryObjective, ") smaller than largest stimulus (", maxValue , ")!")
    return(10^5+abs(maxValue-rightBoundaryObjective)) # Return a large value for convenience for optim
  }
  scaling <- (maxValue - minValue)
  leftBoundary <- minValue - scaling * exp(leftBoundaryExpansion) 
  rightBoundary <- maxValue + scaling * exp(rightBoundaryExpansion) 
  
 
    
 
  if(mode=="Objective Log Likelihood"){
    
    predictions <- stimuli %>% multiCycle(c(leftBoundary, rightBoundary)) %>%  
      psiLogOdds() %>% vanillaBayes(kappa=kappa
                                    , tauStimuli=tauStimuli
                                    , tauCategory= tauCategory
                                    , responses=NULL
                                    , mode="prediction"
      ) %>% psiLogOddsInverse() %>%  multiCycleInverse(c(leftBoundary, rightBoundary)) 
    #plot(predictions)
    0-sum(log(dnorm(predictions-responses, sd=1/(tauStimuli+tauCategory))))
  } else {
     stimuli %>% multiCycle(c(leftBoundary, rightBoundary)) %>%  
      psiLogOdds() %>% vanillaBayes(kappa=kappa
                                    , tauStimuli=tauStimuli
                                    , tauCategory= tauCategory
                                    , responses=(responses  %>% multiCycle(c(leftBoundary, rightBoundary)) %>% psiLogOdds)
                                    , mode=mode
      ) %>% psiLogOddsInverse() %>%  multiCycleInverse(c(leftBoundary, rightBoundary)) 
  }
}


#' fitWarpedBayesModel
#' 
#' A convenience function that packages several commonly popular moves that let a model do optimizaiotn. 
#' @param model a model that has the general layout of the "bayesian..." models included in this package.
#' @param stimuli a vector of stimuli, in whatever raw format you like.
#' @param responses  a vector of stimuli, in whatever raw format you like.
#' @param initialPars an initial set of any parameter values you expect optim to optimize over
#' @param fixedPars a fixed set of any parameter values you do not expect optim to optimize over
#' @param control Passed directly into optim's control parameter
#' @param fit If set to false, simply returns the function evaluation over the initial parameters. This is mostly useful for debugging
#' @param optimizing What criterion should be optimized?  Current valid values are "Objective RMSE" and "Subjective Log Likelihood". Subjective log likelihood is 
#' currently considered to be confusing, and is not recommended for the naive user.
#' @return A tibble that contains one row for each stimulus/response pair, and includes 
#' several columns (see details for details)
#' @details The returned avalue includes two columns that are different on each line--the meanExpecation of the
#' fitted model, and one simulation sampled from the final parameters.  Several more columns pass through the results of hte
#' fit (value, counts, and convergence).  Finally, one column will be made per parameter. 
#' This format is a convenient one if you plan to attach your model fits (and predictions) to a stimulus tibble.
#' @seealso bayesianHuttenlocherSpatialMemory
#' @export
#' @examples
#' a <- fitWarpedBayesModel(bayesianGonzalezWu, 
#'                          fakeStims, 
#'                          fakeData, 
#'                          initialPars = c(kappa=1, tauStimuli=100, tauCategory=10))

fitWarpedBayesModel <- function(model, stimuli, responses
                , initialPars=NULL
                , fixedPars =NULL
                , control=list(maxit=5000, reltol = 10e-200)
                , fit=TRUE
                , optimizing="Objective Log Likelihood"
){
  # Safety check parameters
  if(!(optimizing %in% c("Objective Log Likelihood", "Subjective Log Likelihood"))){
   if(is.character(optimizing)){
     stop("I don't know how to optimize ", optimizing, ". Please give me either Objective RMSE or Subjective Log Likelihood")
   } else {
     stop("Invalid value for optimizing. Please give me either Objective RMSE or Subjective Log Likelihood")
   }
  }

    fitFunction <- function(pars){
      do.call(model, append(append(append(list(stimuli=stimuli), pars), fixedPars), list(responses=responses, mode=optimizing)))
  }
  if(fit){
    result <- stats::optim(initialPars, fitFunction, control=control, method=c("Nelder-Mead") )
    #print(result)
    simulation <- do.call(model, append(append(list(stimuli=stimuli), result$par), list(mode="simulation")))
    meanExpectation <- do.call(model, append(append(list(stimuli=stimuli), result$par), list(mode="prediction")))
  
    a <- tibble::tibble(
      stimulus = stimuli
      , response = responses
      , meanExpectation = meanExpectation
      , simulation = simulation
      , value = result$value
      , counts=result$counts[1]
      , convergence = result$convergence
    )
    for(i in 1:length(result$par)){
      a[names(result$par[i])] <- result$par[i]
    }
    a
  } else { # just for debugging
    result <-fitFunction(initialPars)
    result
  }
}







