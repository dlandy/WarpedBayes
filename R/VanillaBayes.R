
#' vanillaBayes
#' 
#' Performs vanilla bayes (single category) on a set of stimuli
#' @param stimuli a vector of stimuli, between -inf and inf
#' @param kappa The location of the category
#' @param tauStimuli The precision of the stimulus traces: may be a single number or a vector
#' @param tauCategory The precision of the category distribution
#' @param responses an optional vector of responses. If responses are given, the return value is the probability of these responses
#' @return A vector containing mean stimulus locations
#' @details If the kappa, tauStimuli, and tauCategory items are all more than length 1,
#'  and  are length 2 less than the number of bins, then we pad them by negative and positive infinity.
#' @keywords bayesianInference
#' @seealso vanillaBayes
#' @export
#' @examples
#' (0:1000/1000) %>% psiLogOdds() %>% vanillaBayes(kappa=5) %>% psiLogOddsInverse()  # Implements Gonzales & Wu, 1996
#' 1:1000 %>% psiLog() %>% vanillaBayes() %>% psiLogInverse()  # Implements Stevens Power Law
#' plot(-99:100/100, (-99:100/100) %>% multiCycle(references= c(-10, 0, 10)) %>% psiLogOdds() %>% vanillaBayes(kappa=c(-1, 1), tauStimuli=10) %>% psiLogOddsInverse() %>% multiCycleInverse(references=c(-10, 0, 10))-(-99:100/100), ylab="bias", xlab="stimulus");abline(0,0)
vanillaBayes <- function (stimuli, ...) {
  UseMethod("vanillaBayes", stimuli)
}

#' @export
vanillaBayes.list <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1, responses="none"){
  if(length(kappa) == length(stimuli)-2       && length(kappa)>1      ){kappa       <- c(Inf, kappa, Inf)}
  if(length(tauStimuli) == length(stimuli)-2  && length(tauStimuli)>1 ){tauStimuli  <- c(Inf, tauStimuli, Inf)}
  if(length(tauCategory) == length(stimuli)-2 && length(tauCategory)>1){tauCategory <- c(Inf, tauCategory, Inf)}
  mapply(vanillaBayes, stimuli, kappa, tauStimuli, tauCategory, responses, SIMPLIFY=FALSE)
}


#' @export
vanillaBayes.numeric <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1, responses="none") {
  predictions = vanillaBayes.predictions(stimuli, kappa, tauStimuli, tauCategory)
  tauIntegration = tauStimuli + tauCategory
  if(length(responses)==1 && responses=="none"){
    predictions
  } else if(length(responses)==1 && responses=="simulation"){
      rnorm(length(predictions), mean=predictions, sd=1/sqrt(tauIntegration))
  } else {
    if(tauStimuli <= 0 | tauCategory <=0){return(999999)} # large value if tau's go negative
    result <- 0-sum(log(dnorm(predictions-responses, sd=1/sqrt(tauIntegration)))) # Bad Normal Assumption
    class(result) <- append("logLikelihoodOfResponses", class(result))
    result
  }
}

#' @export
vanillaBayes.predictions <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1){
  tauIntegration = tauStimuli + tauCategory
  beta <- tauStimuli/tauIntegration
  beta*stimuli + (1-beta)*kappa
}

