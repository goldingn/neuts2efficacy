#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param x
#' @param c50
#' @param mean_log10_neut
#' @param sd_log10_neut
#' @param log_k
#' @return
#' @author Nick Golding
#' @export
# define the logit-normal density to integrate
logit_normal_density <- function(x, c50, mean_log10_neut, sd_log10_neut, log_k) {
  prob <- prob_avoid_outcome(log10_neut = x, log_k = log_k, c50 = c50)
  dens <- log10_neut_density(x, mean_log10_neut, sd_log10_neut)
  prob * dens
}
