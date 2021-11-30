#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ve
#' @param optimal_log10_neut
#' @return
#' @author Nick Golding
#' @export
# given matching vectors of neut titres and VEs for the same outcome, find
# the c50 that best fits (by least squares)
find_c50 <- function(ves, log10_neuts, c50_interval = c(0.01, 10)) {

  # compute the sum of squares between the expected and observed VEs, for this c50
  optimiser_fun <- function(c50) {

    expected_ves <- vapply(
      log10_neuts,
      FUN = ve_from_mean_log10_neut,
      sd_log10_neut = sd_log10_neut_titres,
      log_k = log_k,
      c50 = c50,
      FUN.VALUE = numeric(1)
    )

    # sum of squares
    sum((ves - expected_ves)^2)

  }

  result <- optimise(optimiser_fun, interval = log10(c50_interval))
  result$minimum

}
