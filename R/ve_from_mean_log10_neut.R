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
  mean_log10_neut_vec,
  sd_log10_neut,
  log_k,
  c50_vec,
  method = c("adaptive", "gaussian"),
  lower = -10,
  upper = 10
) {

  # choose the method and dispatch to the appropriate integration function
  method <- match.arg(method)

  integrator <- switch(
    method,
    adaptive = adaptive_ve_integrator,
    gaussian = gaussian_ve_integrator
  )

  integrals <- integrator(
    c50_vec = c50_vec,
    mean_log10_neut_vec = mean_log10_neut_vec,
    sd_log10_neut = sd_log10_neut,
    log_k = log_k,
    lower = lower,
    upper = upper
  )

  integrals

}
