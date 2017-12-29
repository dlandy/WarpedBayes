
#' psiLinear
#' 
#' Creates mappings from -inf:+inf -> -inf:+inf, with a linear transformation mechanism
#' @param stimuli a vector of stimuli, between -inf and inf
#' @param shift A scalar by which to (positively) offset the scaled values
#' @param scaling A scalar multiplier to the (unshifted) values
#' @return A vector containing warped stimuli
#' @keywords psi perceptual tranformations
#' @seealso psiIdentity, psiLogOdds, psiPrelec, psiLog, psiLinearInverse
#' @export
#' @examples
#' psiLinear(c(0.1, 0.2, 0.3))
#' psiLinear(-10:10)
#' psiLinear(-10:10, shift=10, scaling=2)
#' -10:10 %>% psiLinear(shift=2, scaling=0.5)

