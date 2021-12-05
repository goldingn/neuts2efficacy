#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param neut_model
#' @param draws
#' @return
#' @author Nick Golding
#' @export
# get joint posterior samples of Omicron parameters
sim_omicron_params <- function(neut_model, draws, n = 1e4) {

  R0 <- 6 * neut_model$model_objects$R0_ratio
  titre_fold <- -neut_model$model_objects$omicron_log10_neut_fold

  # get VEs for omicron and delta
  omicron_overall_ve <- 1 - neut_model$model_objects$omicron_vaccine_reduction
  delta_overall_ve <- 1 - neut_model$model_objects$delta_vaccine_reduction

  # % reduction in overall VE against transmission
  immune_evasion <- 1 - omicron_overall_ve / delta_overall_ve


  immune_evasion <- neut_model$model_objects$omicron_vaccine_reduction /
    neut_model$model_objects$delta_vaccine_reduction

  reff_ratio <- neut_model$model_objects$expected_reff_ratio

  sims <- calculate(
    R0, titre_fold, immune_evasion, reff_ratio,
    values = draws,
    nsim = n
  )

  tibble(
    R0 = sims$R0[, 1, 1],
    titre_fold = sims$titre_fold[, 1, 1],
    immune_evasion = sims$immune_evasion[, 1, 1],
    reff_ratio = sims$reff_ratio[, 1, 1]
  )

}
