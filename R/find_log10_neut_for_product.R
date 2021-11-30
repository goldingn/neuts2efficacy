#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ve_symptomatic
#' @param c50_symptoms
#' @return
#' @author Nick Golding
#' @export
# find the mean neutralising antibody titre induced by a product from its VE against
# symptomatic disease (probably the most well understood effect and one that Deb
# has already estimated the c50 for)
find_log10_neut_for_product <- function(ve_symptomatic, c50_symptoms = -0.6966127, log10_neut_interval = log10(c(0.001, 20))) {

  optimiser_fun <- function(log10_neut) {
    expected_ve <- ve_from_mean_log10_neut(
      mean_log10_neut = log10_neut,
      sd_log10_neut = sd_log10_neut_titres,
      log_k = log_k,
      c50 = c50_symptoms
    )
    observed_ve <- ve_symptomatic
    diff <- (observed_ve - expected_ve) ^ 2
    diff
  }

  result <- optimise(f = optimiser_fun, interval = log10_neut_interval)
  log10_neut <- result$minimum
  log10_neut

}
