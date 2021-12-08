#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param neut_model
#' @return
#' @author Nick Golding
#' @export
# get the degree of immune evasion (among the vaccinated) from the neut model
immune_evasion <- function(neut_model) {

  za_baseline_immunity_log10_neut_fold <- neut_model$model_objects$za_baseline_immunity_log10_neut_fold
  omicron_log10_neut_fold <- neut_model$model_objects$omicron_log10_neut_fold

  omicron_ves_peak <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = rep(za_baseline_immunity_log10_neut_fold + omicron_log10_neut_fold, 5),
    sd_log10_neut = neut_model$model_objects$sd_log10_neut_titres,
    log_k = neut_model$model_objects$log_k,
    c50_vec = neut_model$model_objects$c50s,
    method =  "gaussian"
  )

  delta_ves_peak <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = rep(za_baseline_immunity_log10_neut_fold, 5),
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
  1 - omicron_overall_ve / delta_overall_ve

}

