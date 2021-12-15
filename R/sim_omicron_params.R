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
sim_omicron_params <- function(neut_model, draws,
                               n = 1e4,
                               baseline_immunity = "Pfizer_dose_2") {

  # recompute VEs against transmission for Omicron and Delta, to remove the
  # effect of waning
  za_distancing_effect <- neut_model$model_objects$za_distancing_effect
  omicron_log10_neut_fold <- neut_model$model_objects$omicron_log10_neut_fold

  evasion <- immune_evasion(neut_model, baseline_immunity = baseline_immunity)

  # other ratios
  R0 <- 6 * neut_model$model_objects$R0_ratio
  titre_fold <- -neut_model$model_objects$omicron_log10_neut_fold

  reff_ratio <- neut_model$model_objects$expected_reff_ratio

  R0_ratio <- neut_model$model_objects$R0_ratio

  sims <- calculate(
    R0, titre_fold, evasion, reff_ratio, R0_ratio, za_distancing_effect,
    values = draws,
    nsim = n
  )

  tibble(
    R0 = sims$R0[, 1, 1],
    titre_fold = sims$titre_fold[, 1, 1],
    immune_evasion = sims$evasion[, 1, 1],
    reff_ratio = sims$reff_ratio[, 1, 1],
    R0_ratio = sims$R0_ratio[, 1, 1],
    za_distancing_effect = sims$za_distancing_effect[, 1, 1]
  )

}
