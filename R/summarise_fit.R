#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param draws
#' @return
#' @author Nick Golding
#' @export
summarise_fit <- function(draws) {

  rhats <- coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)
  neffs <- coda::effectiveSize(draws)
  summaries <- summary(draws)

  tibble(
    parameter = names(neffs),
    mean = summaries$statistics[, "Mean"],
    sd = summaries$statistics[, "SD"],
    lower_95_ci = summaries$quantiles[, "2.5%"],
    upper_95_ci = summaries$quantiles[, "97.5%"],
    rhat = rhats$psrf[, "Point est."],
    rhat_upper = rhats$psrf[, "Upper C.I."],
    n_eff = neffs,
  )

}
