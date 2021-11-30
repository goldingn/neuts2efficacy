#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param log10_neut
#' @param log_k
#' @param c50
#' @return
#' @author Nick Golding
#' @export
# logistic function or the individual's probability of avoiding the outcome
# (being infected, passing on infection, symptomatic disease, severe disease,
# death) given their log10 neut titre, the C50 for that outcome, and the
# modelled slop of the relationship
prob_avoid_outcome <- function(log10_neut, log_k, c50) {
  1 / (1 + exp(-exp(log_k) * (log10_neut - c50)))
}
