#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
get_far_sar_coefs <- function() {

  # pairs of household FAR and SAR values for households with one source
  # infection and m susceptibles, from reading off Fig 1 in Sharker et al.
  # https://doi.org/10.1371/journal.pcbi.1008601
  far_sar <- tribble(
    ~m, ~far, ~sar,
    2, 0.11, 0.1,
    2, 0.37, 0.3,
    4, 0.13, 0.1,
    4, 0.52, 0.3,
  )

  # m = 3 is about right for Lyngse et al data, based on distributionof
  # susceptibles for both variants, that's not plotted so linearly interpolate m
  # = 2 and m = 4 by fitting to both
  m <- lm(far ~ sar, data = far_sar)

  m$coefficients

}
