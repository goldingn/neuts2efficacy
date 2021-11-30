#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param mean_log10_neut
#' @param sd_log10_neut
#' @param log_k
#' @param c50
#' @return
#' @author Nick Golding
#' @export
# given a neut population mean neut titre, standard deviation in the population
# neut titres, slop parameter and C50 for that outcome, compute the
# population-level VE by integrating over population neut titres
ve_from_mean_log10_neut <- function(
  mean_log10_neut,
  sd_log10_neut,
  log_k,
  c50
) {

  integrate_fun <- function(x) {
    prob <- prob_avoid_outcome(log10_neut = x, log_k = log_k, c50 = c50)
    weight <- dnorm(x, mean_log10_neut, sd_log10_neut)
    prob * weight
  }

  integral <- integrate(
    f = integrate_fun,
    lower = mean_log10_neut - 5 * sd_log10_neut,
    upper = mean_log10_neut + 5 * sd_log10_neut
  )

  integral$value

}
