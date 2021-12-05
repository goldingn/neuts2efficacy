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

  # recompute VEs against transmission for Omicron and Delta, to remove the
  # effect of waning

  omicron_log10_neut_fold <- neut_model$model_objects$omicron_log10_neut_fold

  omicron_ves_peak <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = rep(omicron_log10_neut_fold, 5),
    sd_log10_neut = neut_model$model_objects$sd_log10_neut_titres,
    log_k = neut_model$model_objects$log_k,
    c50_vec = neut_model$model_objects$c50s,
    method =  "gaussian"
  )

  delta_ves_peak <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = as.matrix(rep(0, 5)),
    sd_log10_neut = neut_model$model_objects$sd_log10_neut_titres,
    log_k = neut_model$model_objects$log_k,
    c50_vec = neut_model$model_objects$c50s,
    method =  "gaussian"
  )

  acquisition_idx <- match("acquisition", neut_model$lookups$outcome)
  transmission_idx <- match("transmission", neut_model$lookups$outcome)

  omicron_vaccine_reduction <- prod(1 - omicron_ves_peak[c(acquisition_idx, transmission_idx)])
  delta_vaccine_reduction <- prod(1 - delta_ves_peak[c(acquisition_idx, transmission_idx)])

  # VEs against overall transmission for omicron and delta
  omicron_overall_ve <- 1 - omicron_vaccine_reduction
  delta_overall_ve <- 1 - delta_vaccine_reduction

  # % reduction in overall VE against transmission
  immune_evasion <- 1 - omicron_overall_ve / delta_overall_ve

  # other ratios
  R0 <- 6 * neut_model$model_objects$R0_ratio
  titre_fold <- -neut_model$model_objects$omicron_log10_neut_fold

  reff_ratio <- neut_model$model_objects$expected_reff_ratio

  R0_ratio <- neut_model$model_objects$R0_ratio

  sims <- calculate(
    R0, titre_fold, immune_evasion, reff_ratio, R0_ratio,
    values = draws,
    nsim = n
  )

  tibble(
    R0 = sims$R0[, 1, 1],
    titre_fold = sims$titre_fold[, 1, 1],
    immune_evasion = sims$immune_evasion[, 1, 1],
    reff_ratio = sims$reff_ratio[, 1, 1],
    R0_ratio = sims$R0_ratio[, 1, 1]
  )

}
