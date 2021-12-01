#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param c50_vec
#' @param mean_log10_neut_vec
#' @param sd_log10_neut
#' @param log_k
#' @param lower
#' @param upper
#' @return
#' @author Nick Golding
#' @export
adaptive_ve_integrator <- function(
  c50_vec,
  mean_log10_neut_vec,
  sd_log10_neut,
  log_k,
  lower,
  upper
) {


  if (
    inherits(c50_vec, "greta_array") |
    inherits(mean_log10_neut_vec, "greta_array") |
    inherits(sd_log10_neut, "greta_array") |
    inherits(log_k, "greta_array")
  ) {
    stop ("adaptive integration can not be used with greta models")
  }

  integrate_once <- function(c50, mean_log10_neut) {
    integral <- integrate(
      f = logit_normal_density,
      c50 = c50,
      mean_log10_neut = mean_log10_neut,
      sd_log10_neut = sd_log10_neut,
      log_k = log_k,
      lower = lower,
      upper = upper
    )$value
  }

  logit_normal_density(0, c50_vec[1], mean_log10_neut = mean_log10_neut_vec[1], sd_log10_neut = sd_log10_neut, log_k = log_k)

  integrals <- mapply(
    integrate_once,
    c50_vec,
    mean_log10_neut_vec
  )

  integrals

}
