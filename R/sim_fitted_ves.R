#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param neut_model
#' @return
#' @author Nick Golding
#' @export
sim_fitted_ves <- function(neut_model, draws) {

  ve_expected <- neut_model$model_objects$ve_expected

  sims <- calculate(
    ve_expected,
    values = draws,
    nsim = 1000
  )

  sims[[1]][, , 1]
}
