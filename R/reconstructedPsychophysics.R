# Reconstructed classical functions, using Bayesian parameterizations

robustLog <- function(x, smallValue=10^-322){
  x[x<=smallValue] <- smallValue
  log(x)
} 


negSumLogs <- function(probabilities){
  return(0-sum(robustLog(probabilities)))
}







#' uniformGuessingModel
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param responses an optional vector of responses
#' @param responseGrid An optional vector of possible responses
#' @param mode What aspect should the function calculate? Legel choices include "prediction", "simulation", and "subjectiveLogLikelihood"
#' If responses are given, the return value is the logLikelihood of the responses given the parameters
#' @return A vector the transformed stimuli
#' @seealso bayesianHuttenlocherSpatialMemory
#' @export
#' @examples
#' uniformGuessingModel(1:100/100)

uniformGuessingModel <- function(stimuli
                                 , responses=NULL
                                 , responseGrid = NULL
                                 , mode = "prediction"
){
  nStim <- length(stimuli)
  if(is.null(responseGrid)){

    warning("You didn't give me a grid of responses. I made one up, but that's probably not very safe.")
    responseGrid <- unlist(ifelse(is.null(responses)
                           , ifelse(is.null(stimuli), list(c(0)), list(sort(unique(stimuli))))
                           , list(sort(unique(responses)))
    )
                           
    )
  }
  if(mode=="prediction"){
    return(rep(median(responseGrid), nStim))
  } else if(mode=="simulation"){
    return(sample(responseGrid, nStim, replace=TRUE))
  } else if(mode=="subjectiveLogLikelihood"){
    return(rep(1/length(responseGrid), nStim))
  }
}

#' guessingModel
#' 
#' This function takes a model, and returns a version of that model that incorporates uniform guessing at some 'guessRate'.
#' It is very important to give guessing models a responseGrid parameter (otherwise, they don't know how to distribute guess probabilities!).
#' Otherwise, these models can be used essentially in the same ways as normal models
#' 
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param responses an optional vector of responses
#' @param responseGrid An optional vector of possible responses
#' @param mode What aspect should the function calculate? Legel choices include "prediction", "simulation", and "subjectiveLogLikelihood"
#' If responses are given, the return value is the logLikelihood of the responses given the parameters
#' @param ... All parameters that should be passed into the main model
#' @return A vector the transformed stimuli
#' @seealso bayesianHuttenlocherSpatialMemory
#' @export
#' @examples
#' (guessingModel(bayesianGonzalezWu))(1:9/10, mode="prediction", guessRate=0.5,  responseGrid=1:9/10)
#' a <- with(boundedProportionWithGuessingSimulatedData,
#' fitWarpedBayesModel(guessingModel(bayesianGonzalezWu), 
#'                     stimulus
#'                     ,  response
#'                     , initialPars = c(guessRate = 0, kappa=1, tauStimuli=10, tauCategory=10, leftBoundaryExpansion=-1, rightBoundaryExpansion=-1)
#'                     , responseGrid = sort(unique(response))
#'                     , fixedPars=c() ))
guessingModel <- function(model){
  function(
                           stimuli=NULL
                          , responses=NULL
                          , responseGrid = sort(unique(responses))
                          , mode = "prediction"
                          , guessRate=0
                          , ...){
    if(mode=="simulation"){
      guess <- uniformGuessingModel(stimuli=stimuli
                                    , responses=responses
                                    , responseGrid=responseGrid
                                    , mode=mode)
      nonGuess <- model(stimuli=stimuli
                        , responses=responses
                        , responseGrid=responseGrid
                        , mode=mode, ...)
      pick <- rbinom(length(stimuli), 1,guessRate)
      return(pick*guess+(1-pick)*nonGuess)
      
      
    } else if(mode=="subjectiveLogLikelihood" | mode=="prediction" | TRUE){ #Secretly a catch-all
    guessRate*uniformGuessingModel(stimuli=stimuli
                                   , responses=responses
                                   , responseGrid=responseGrid
                                   , mode=mode) + 
    (1-guessRate)*model(stimuli=stimuli
                        , responses=responses
                        , responseGrid=responseGrid
                        , mode=mode, ...)
    
  }
    
}
}



#' fitWarpedBayesModel
#' 
#' A convenience function that packages several commonly popular moves that let a model do optimizaiotn. 
#' @param model a model that has the general layout of the "bayesian..." models included in this package.
#' @param stimuli a vector of stimuli, in whatever raw format you like.
#' @param responses  a vector of stimuli, in whatever raw format you like.
#' @param responseGrid a vector of possible responses, for discretizing the normal distribution
#' @param initialPars an initial set of any parameter values you expect optim to optimize over
#' @param fixedPars a fixed set of any parameter values you do not expect optim to optimize over
#' @param control Passed directly into optim's control parameter
#' @param fit If set to false, simply returns the function evaluation over the initial parameters. This is mostly useful for debugging
#' @param optimizing What criterion should be optimized?  Current valid values are "Objective RMSE" and "subjectiveLogLikelihood". subjectiveLogLikelihood is 
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
#' fitWarpedBayesModel(bayesianSpatialMemoryLandyCrawfordCorbin2017, 
#'                     stimulus
#'                     ,  response
#'                     , initialPars = c(kappa=1, tauStimuli=10, tauCategory=10, leftBoundaryExpansion=-1, rightBoundaryExpansion=-1)
#'                     , responseGrid = sort(unique(response))
#'                     , fixedPars=c() ))

fitWarpedBayesModel <- function(model, stimuli, responses
                , responseGrid = NULL
                , initialPars=NULL
                , fixedPars =NULL
                , control=list(maxit=5000, reltol = 10e-200)
                , fit=TRUE
                , optimizing="subjectiveLogLikelihood"
){
  # Safety check parameters
  if(!(optimizing %in% c("subjectiveLogLikelihood"))){
   if(is.character(optimizing)){
     stop("I don't know how to optimize ", optimizing, ". Please give me subjectiveLogLikelihood")
   } else {
     stop("Invalid value for optimizing. Please give me either Objective RMSE or subjectiveLogLikelihood")
   }
  }
  if(is.null(responseGrid)){
    warning("You didn't give me a grid of responses. I made one up, but that's probably not very safe.")
    responseGrid = sort(unique(responses))
  }
    fitFunction <- function(pars){
      negSumLogs(do.call(model, append(append(append(list(stimuli=stimuli), pars), fixedPars), 
                                       list(responses=responses, mode=optimizing, responseGrid=responseGrid))))

      
  }
  if(fit){
    result <- stats::optim(initialPars, fitFunction, control=control, method=c("Nelder-Mead") )
    simulation <- do.call(model, append(append(append(list(stimuli=stimuli), result$par), fixedPars), 
                                        list(mode="simulation", responseGrid=responseGrid)))
    meanExpectation <- do.call(model, append(append(append(list(stimuli=stimuli), result$par), fixedPars), list(mode="prediction", responseGrid=responseGrid)))
    #print(result)
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









#' bayesianGonzalezWu
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param responses an optional vector of responses
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
#' @param minValue The lowest range value (in objective units).  Should be at least as small as the smallest stimulus and response.
#' @param maxValue The highest range value  (in objective units).  Should be at least as large as the largest stimulus and response.
#' @param smallValue a small amount, by which to default boundary expansion
#' @param mode What aspect should the function calculate? Legel choices include "prediction", "simulation", and "subjectiveLogLikelihood"
#' @param responseGrid an optional vector of response structured Responses
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
                               , responses=NULL
                               , tauStimuli = 1
                               , tauCategory = 1
                               , leftBoundaryObjective  = minValue-smallValue
                               , rightBoundaryObjective = maxValue + smallValue
                               , rightBoundaryExpansion = NULL
                               , leftBoundaryExpansion  = NULL
                               , minValue = min(c(stimuli, responses), na.rm=T)
                               , maxValue = max(c(stimuli, responses), na.rm=T)
                               , smallValue = 10^-10
                               , mode = "prediction"
                               , responseGrid = NULL
                               ){
  if(is.null(responseGrid)){
    if(mode %in% c("subjectiveLogLikelihood")){ 
      warning("You didn't give me a grid of responses. I made one up, but that's probably not very safe.")
    }
    responseGrid <- unlist( ifelse(is.null(responses), list(c(0)), list(sort(unique(responses)))))
  }
  leftBoundary  <- ifelse(is.null(leftBoundaryExpansion ), leftBoundaryObjective,   minValue -  exp(leftBoundaryExpansion ) )
  rightBoundary <- ifelse(is.null(rightBoundaryExpansion), rightBoundaryObjective,  maxValue +  exp(rightBoundaryExpansion) )
  
  perception <- function(vals) { vals %>% uniCycle(c(leftBoundary, rightBoundary)) %>%  psiLogOdds()}
  
  
  a <- perception(stimuli) %>% 
    vanillaBayes(kappa=kappa
                 , tauStimuli=tauStimuli
                 , tauCategory= tauCategory
                 , responses=perception(responses)
                 , mode=mode
                 , responseGrid = perception(responseGrid)
    ) %>% 
    psiLogOddsInverse() %>%  
    uniCycleInverse(c(leftBoundary, rightBoundary)) 
  
  a
}





#' bayesianStevensPowerLaw
#' @param stimuli a vector of stimuli, between 0 and inf
#' @param kappa The location of the category
#' @param kappaObjective An alternative specification giving kappa in objective units
#' @param tauStimuli The precision of the stimulus traces: may be a single number or a vector
#' @param tauCategory The precision of the category distribution
#' @param responses an optional vector of responses.
#' @param responseGrid an optional vector of response structured Responses
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
                                    , kappa=psiLog(kappaObjective), tauStimuli=1, tauCategory=1, responses=NULL
                                    , mode ="prediction"
                                    , responseGrid = NULL
){
  if(is.null(responseGrid) & mode %in% c("subjectiveLogLikelihood")){
    warning("You didn't give me a grid of responses. I made one up, but that's probably not very safe.")
    responseGrid <- unlist( ifelse(is.null(responses), list(c(0)), list(sort(unique(responses)))))
  }

    stimuli %>% psiLog %>% vanillaBayes(kappa=kappa
                                        , tauStimuli=tauStimuli
                                        , tauCategory=tauCategory
                                        , responses=psiLog(responses)
                                        , responseGrid = psiLog(responseGrid)
                                        ) %>% psiLogInverse()
}


#' bayesianSpatialMemoryHuttenlocher
#' 
#' A simple package that assembles the symmetric model used by Huttenlocher and colleagues to analyze spatial 
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
#' @param responseGrid an optional vector of response structured Responses
#' @return A vector the transformed stimuli, or the logLikelihood of them.
#' @details This package 
#' @seealso psiIdentity, multiCycleInverse
#' @export
#' @examples
#' bayesianSpatialMemoryHuttenlocher(-99:100/100)
#' bayesianSpatialMemoryHuttenlocher(-99:100/100, kappa=1, tauStimuli=2)
#' bayesianSpatialMemoryHuttenlocher(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianSpatialMemoryHuttenlocher <- function(stimuli
                                              , responses=NULL
                                              , kappaObjective = 0.5
                                              , kappa=psiLogOdds(kappaObjective)
                                              , tauStimuli=1
                                              , tauCategory=1
                                              , boundary = 1
                                              , leftBoundary = -1*boundary
                                              , rightBoundary = boundary
                                              , center = 0
                                              , mode="prediction"
                                              , responseGrid = NULL
                                              ){
  if(is.null(responseGrid) & mode %in% c("subjectiveLogLikelihood")){
    warning("You didn't give me a grid of responses. I made one up, but that's probably not very safe.")
    responseGrid <- unlist( ifelse(is.null(responses), list(c(0)), list(sort(unique(responses)))))
  }
  refs <- c(leftBoundary, center, rightBoundary)
  stimuli %>% multiCycle(references = refs) %>%
    vanillaBayes(kappa=kappa
                 , tauStimuli=tauStimuli
                 , tauCategory=tauCategory
                 , responses   = multiCycle(responses,    references = refs)
                 , responseGrid= multiCycle(responseGrid, references = refs)
                 , mode=mode
                 )  %>% 
    multiCycleInverse(references = refs)
}




#' bayesianSpatialMemoryLandyCrawfordCorbin2017
#' 
#' A simple package that assembles the symmetric model used by Landy, Crawford, and Corbin, 2017 to analyze spatial 
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
#' @param responseGrid an optional vector of response structured Responses
#' @return A vector the transformed stimuli, or the logLikelihood of them.
#' @seealso psiIdentity, multiCycleInverse
#' @export
#' @examples
#' bayesianSpatialMemoryLandyCrawfordCorbin2017(-99:100/100)
#' bayesianSpatialMemoryLandyCrawfordCorbin2017(-99:100/100, kappa=1, tauStimuli=2)
#' bayesianSpatialMemoryLandyCrawfordCorbin2017(1:100, kappa=1, tauStimuli=2, responses=2*(1:100)^.9)
bayesianSpatialMemoryLandyCrawfordCorbin2017 <- function(stimuli
                                                         , kappa = psiLogOdds(multiCycle(kappaObjective, references=c(leftBoundaryObjective, rightBoundaryObjective)))
                                                         , kappaObjective = 0.5
                                                         , tauStimuli = 1
                                                         , tauCategory = 1
                                                         , leftBoundaryObjective  = minValue-smallValue
                                                         , rightBoundaryObjective = maxValue+smallValue
                                                         , rightBoundaryExpansion = NULL
                                                         , leftBoundaryExpansion  = NULL
                                                         , minValue = min(c(stimuli), na.rm=T)
                                                         , maxValue = max(c(stimuli), na.rm=T)
                                                         , smallValue = 10^-10
                                                         , responses=NULL
                                                         , center = 0
                                                         , mode = "prediction"
                                                         , responseGrid = NULL
                                                         ){
  if(is.null(responseGrid) & mode %in% c("subjectiveLogLikelihood")){
    warning("You didn't give me a grid of responses. I made one up, but that's probably not very safe.")
    responseGrid <- unlist( ifelse(is.null(responses), list(c(0)), list(sort(unique(responses)))))
  }
  leftBoundary  <- ifelse(is.null(leftBoundaryExpansion ), leftBoundaryObjective,   minValue -  exp(leftBoundaryExpansion ) )
  rightBoundary <- ifelse(is.null(rightBoundaryExpansion), rightBoundaryObjective,  maxValue +  exp(rightBoundaryExpansion) )
  refs <- c(leftBoundary, center, rightBoundary)
  kappas <- c(kappa, 1-kappa)
  stimuli %>% multiCycle(references = refs) %>%
      psiLogOdds() %>% 
      vanillaBayes(kappa=kappas
                   , tauStimuli=tauStimuli
                   , tauCategory=tauCategory
                   , responses    = multiCycle(responses,    references = refs) %>% psiLogOdds()
                   , responseGrid = multiCycle(responseGrid, references = refs) %>% psiLogOdds()
                   , mode=mode
      )  %>% 
      psiLogOddsInverse() %>% 
      multiCycleInverse(references = refs)
  
}


