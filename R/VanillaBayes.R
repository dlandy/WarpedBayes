
#' vanillaBayes
#' 
#' Performs vanilla bayes (single category) on a set of stimuli
#' @param stimuli a vector of stimuli, between -inf and inf
#' @param kappa The location of the category
#' @param tauStimuli The precision of the stimulus traces: may be a single number or a vector
#' @param tauCategory The precision of the category distribution
#' @param responses an optional vector of responses. If responses are given, the return value is the 
#' @return A vector containing mean stimulus locations
#' @keywords bayesianInference
#' @seealso vanillaBayes
#' @export
#' @examples
#' (0:1000/1000) %>% psiLogOdds() %>% vanillaBayes(kappa=5) %>% psiLogOddsInverse()  # Implements Gonzales & Wu, 1996
#' 1:1000 %>% psiLog() %>% vanillaBayes() %>% psiLogInverse()  # Implements Stevens Power Law
vanillaBayes <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1, responses="none") {
  predictions = vanillaBayes.predictions(stimuli, kappa, tauStimuli, tauCategory)
  tauIntegration = tauStimuli + tauCategory
  if(responses=="none"){
    predictions
  } else {
    if(tauStimuli <= 0 | tauCategory <=0){return(999999)} # large value if tau's go negative
    0-sum(log(dnorm(predictions-responses, sd=tauIntegration))) # Bad Normal Assumption
  }
}

#' @export
vanillaBayes.predictions <- function(stimuli, kappa=0, tauStimuli=1, tauCategory=1){
  tauIntegration = tauStimuli + tauCategory
  beta <- tauStimuli/tauIntegration
  beta*stimuli + (1-beta)*kappa
}

