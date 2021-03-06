
#' vanillaBayes
#' 
#' Performs vanilla bayes (single category) on a set of stimuli. 
#' @param stimuli a vector of stimuli, between -inf and inf
#' @param kappa The location of the category
#' @param tauStimuli The precision of the stimulus traces: may be a single number or a vector
#' @param tauCategory The precision of the category distribution
#' @param responses an optional vector of responses. Should only be given if mode is "subjectiveLogLikelihood"
#' @param mode What aspect should the function calculate? Legel choices include "prediction", "simulation", and "subjectiveLogLikelihood"
#' @return A vector containing mean predicted stimulus locations, or the log likelihood of the responses given the model
#' @details This function assumes that the data are in a metric 
#' space (-inf, inf), with a single normally distributed generating category (with mean kappa and precision tauCategory). It further assumes a set of stimuli, which are
#' normal distributions with means at the value of stimuli, and precision tauStimuli.  It returns either the mean expected 
#' location of response to the stimuli (given the parameters), or if a set of responses is given, the log likelihood of the 
#' responses given the model.
#' 
#' If the kappa, tauStimuli, and tauCategory items are all more than length 1,
#'  and  are length 2 less than the number of bins, then we pad them by negative and positive infinity.
#' @keywords bayesianInference
#' @seealso vanillaBayes
#' @export
#' @examples
#' (0:1000/1000) %>% vanillaBayes(kappa=5) 
#'     # The  Bayesian normal-normal model typical to many analyses
#' (0:1000/1000) %>% psiLogOdds() %>% 
#'     vanillaBayes(kappa=5) %>% 
#'     psiLogOddsInverse()  #Gonzales & Wu, 1996
#' 1:1000 %>% psiLog() %>% vanillaBayes() %>% psiLogInverse()  # Stevens Power Law
#' plot(-99:100/100, (-99:100/100) %>% multiCycle(references= c(-10, 0, 10)) %>% 
#'     psiLogOdds() %>% vanillaBayes(kappa=c(-1, 1), tauStimuli=10) %>% 
#'     psiLogOddsInverse() %>% 
#'     multiCycleInverse(references=c(-10, 0, 10))-(-99:100/100), 
#'         ylab="bias", xlab="stimulus");abline(0,0)
vanillaBayes <- function (stimuli, kappa=0, tauStimuli=1, tauCategory=1, responses=NULL, mode="prediction", responseGrid = c(0)) {
  UseMethod("vanillaBayes", stimuli)
}

#' @export 
vanillaBayes.list <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1, responses=NULL, mode="prediction", responseGrid = c(0)){
  if(length(kappa) == length(stimuli)-2       && length(kappa)>1      ){kappa       <- c(Inf, kappa, Inf)}
  if(length(tauStimuli) == length(stimuli)-2  && length(tauStimuli)>1 ){tauStimuli  <- c(Inf, tauStimuli, Inf)}
  if(length(tauCategory) == length(stimuli)-2 && length(tauCategory)>1){tauCategory <- c(Inf, tauCategory, Inf)}
  
  if(length(responses)>0){
    result <- mapply(vanillaBayes, stimuli, kappa, tauStimuli, tauCategory, responses, mode=mode, responseGrid = responseGrid, SIMPLIFY=FALSE)
  } else {
    result <- mapply(vanillaBayes, stimuli, kappa, tauStimuli, tauCategory, mode=mode, SIMPLIFY=FALSE)
  }
  # if((length(responses)==0)&&(mode=="prediction" || mode=="simulation" )){
  if((mode=="prediction" || mode=="simulation" )){
    result
  }else{
    result <- (unlist(result))
    
    class(result) <- append("subjectiveLogLikelihood", class(result))
    result
  } 
}

discreteNorm <- 
  function(mean, sd, responseIndex, responseGrid){
    if(length(na.omit(responseGrid[responseIndex+1]))>0){
      
      if(max(is.na(responseGrid[responseIndex+1]))){
        warning("some values of the responseGrid[responseIndex+1] were NaN!", responseGrid[responseIndex+1])
        return(rep(0, length(responseGrid[responseIndex+1])))
        
      }  else if(max(is.na(responseGrid[responseIndex]))){
        warning("some values of the responseGrid[responseIndex] were NaN!", responseGrid[responseIndex])
        return(rep(0, length(responseGrid[responseIndex+1])))
        
      } 
    }
    replaceValues <- is.infinite(mean)
    mean[replaceValues] <- 3e200*sign(mean[replaceValues]) # Sometimes infinite values in the mean throw off pnorm
    a <- pnorm(responseGrid[responseIndex+1]
          , mean=mean
          , sd=sd
          )-
    pnorm(responseGrid[responseIndex]
          , mean=mean
          , sd=sd
          )
    
    a
    
}


#' @export 
vanillaBayes.numeric <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1, responses=NULL, mode="prediction", responseGrid = c(0)) {

  if((!is.null(responses)) && ((mode=="prediction")||(mode=="simulation"))){
    warning("You gave me responses....but you set the mode for ", mode, ".  Did you mean to do this?")
  } 
  vanillaBayesPredictions <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1){
    tauIntegration = tauStimuli + tauCategory
    beta <- tauStimuli/tauIntegration
    beta*stimuli + (1-beta)*kappa
  }
  predictions = vanillaBayesPredictions(stimuli, kappa, tauStimuli, tauCategory)
  tauIntegration = tauStimuli + tauCategory
  if(mode=="prediction"){
    predictions
  } else if(mode=="simulation"){
      rnorm(length(predictions), mean=predictions, sd=1/sqrt(tauIntegration))
  } else if(mode=="subjectiveLogLikelihood") {
    if(tauStimuli <= 0 | tauCategory <=0){# large value if tau's go negative
      result <-rep(0, length(predictions))
    }  else {
      responseIndex <- match(responses, responseGrid)
      responseGrid2 <- (c(-Inf, responseGrid)+c(responseGrid, Inf))/2
    
      #result <- dnorm(predictions-responses, sd=1/sqrt(tauIntegration))# DNORM is not well-normalized here
      result <- discreteNorm(predictions, sd=1/sqrt(tauIntegration), responseIndex, responseGrid2)
     
    }
    class(result) <- append("subjectiveLogLikelihood", class(result))
    result
  }
}




