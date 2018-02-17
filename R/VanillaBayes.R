
#' vanillaBayes
#' 
#' Performs vanilla bayes (single category) on a set of stimuli. 
#' @param stimuli a vector of stimuli, between -inf and inf
#' @param kappa The location of the category
#' @param tauStimuli The precision of the stimulus traces: may be a single number or a vector
#' @param tauCategory The precision of the category distribution
#' @param responses an optional vector of responses. Should only be given if mode is "likelihoodOfResponses"
#' @param mode What aspect should the function calculate? Legel choices include "prediction", "simulation", and "likelihoodOfResponses"
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
vanillaBayes <- function (stimuli, kappa=0, tauStimuli=1, tauCategory=1, responses=NULL, mode="prediction") {
  UseMethod("vanillaBayes", stimuli)
}

#' @export 
vanillaBayes.list <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1, responses=NULL, mode="prediction"){
  if(length(kappa) == length(stimuli)-2       && length(kappa)>1      ){kappa       <- c(Inf, kappa, Inf)}
  if(length(tauStimuli) == length(stimuli)-2  && length(tauStimuli)>1 ){tauStimuli  <- c(Inf, tauStimuli, Inf)}
  if(length(tauCategory) == length(stimuli)-2 && length(tauCategory)>1){tauCategory <- c(Inf, tauCategory, Inf)}
  
  if(length(responses)>0){
    result <- mapply(vanillaBayes, stimuli, kappa, tauStimuli, tauCategory, responses, mode=mode, SIMPLIFY=FALSE)
  } else {
    result <- mapply(vanillaBayes, stimuli, kappa, tauStimuli, tauCategory, mode=mode, SIMPLIFY=FALSE)
  }
  if((length(responses)==0)&&(mode=="prediction" || mode=="simulation" )){
    result
  }else{
    #result <- sum(unlist(result))
    result <- (unlist(result))
    
    class(result) <- append("likelihoodOfResponses", class(result))
    result
  } 
}


#' @export 
vanillaBayes.numeric <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1, responses=NULL, mode="prediction") {

  if((!is.null(responses)) && ((mode=="prediction")||(mode=="simulation"))){
    warning("You gave me responses, so I gave you a log likelihood....but you set the mode for ", mode, ".  Did you mean to do this?")
  } 
  vanillaBayesPredictions <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1){
    tauIntegration = tauStimuli + tauCategory
    beta <- tauStimuli/tauIntegration
    beta*stimuli + (1-beta)*kappa
  }
  robustLog <- function(x, smallValue=10^-300){
    x[x<=smallValue] <- smallValue
    log(x)
  } 
  predictions = vanillaBayesPredictions(stimuli, kappa, tauStimuli, tauCategory)
  tauIntegration = tauStimuli + tauCategory
  if(mode=="prediction"){
    predictions
  } else if(mode=="simulation"){
      rnorm(length(predictions), mean=predictions, sd=1/sqrt(tauIntegration))
  } else if(mode=="likelihoodOfResponses") {
    if(tauStimuli <= 0 | tauCategory <=0){return(999999)} # large value if tau's go negative
    
    
    result <- dnorm(predictions-responses, sd=1/sqrt(tauIntegration))
    class(result) <- append("likelihoodOfResponses", class(result))
    
    result
  }
}




