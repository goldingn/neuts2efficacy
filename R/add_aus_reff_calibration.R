#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param neut_model
#' @return
#' @author Nick Golding
#' @export
add_aus_reff_calibration <- function(neut_model, immunity_strata) {

  # 1. calculate VEs against infection and onward transmission for all combinations
  # of immunity and times since immune event

  # get unwaned neuts for all immunity types

  # copy over code from predict_ves to make this work - or even better, refactor
  # some of it into a separaate function to avoid duplicating it.

  # log10_neuts_list <- list(
  #   neut_model$model_objects$peak_mean_log10_neuts[1],
  #   neut_model$model_objects$peak_mean_log10_neuts[2],
  #   neut_model$model_objects$peak_mean_log10_neuts[3],
  #   neut_model$model_objects$peak_mean_log10_neuts[4],
  #   neut_model$model_objects$peak_mean_log10_neuts[5]
  # )
  #
  # names(log10_neuts_list) <- neut_model$lookups$immunity
  #
  # # pull out the appropriate peak neuts
  # index <- match(immunity_strata$immunity, names(log10_neuts_list))
  # peak_log10_neuts <- log10_neuts_list[index]

  # apply waning to get mean log10 neuts for each stratum
  mean_log10_neuts <- log10_neut_over_time(
    time = immunity_strata$time_since_immunity,
    maximum_log10_neut = peak_log10_neuts,
    decay = neut_model$model_objects$neut_decay
  )

  acquisition_idx <- match("acquisition", neut_model$lookups$outcome)
  transmission_idx <- match("transmission", neut_model$lookups$outcome)

  # calculate all the VEs against infection and onward transmission from these neuts
  all_ves_infection <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = mean_log10_neuts,
    sd_log10_neut = neut_model$model_objects$sd_log10_neut_titres,
    log_k = neut_model$model_objects$log_k,
    c50_vec = neut_model$model_objects$c50s[acquisition_idx],
    method =  "gaussian"
  )

  all_ves_onward <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = mean_log10_neuts,
    sd_log10_neut = neut_model$model_objects$sd_log10_neut_titres,
    log_k = neut_model$model_objects$log_k,
    c50_vec = neut_model$model_objects$c50s[transmission_idx],
    method =  "gaussian"
  )

  # 2. combine these with fractions of the population in each of those bins to get
  # average VEs by age

  # expand these out to infection only, vaccine only, infection and vaccine, etc.
  ves_infection_by_age <- weight_all_ves(all_ves$infection, population_fractions)
  ves_onward_by_age <- weight_all_ves(all_ves$onward, population_fractions)

  # combine with the fractions of the population in that age group that have any
  # immunity, to get multiplliers on infection due to blocking of infection or
  # onward transmission

  # we'll need to expand this out to match the method for combining vaccine and
  # immunity derived effects from the TP/Reff model: compute coverage and VEs of
  # the infected-and-not-vaccinated, the vaccinated-and-not-infected, and the
  # vaccinated-and-infected, and then combine them:


  coverage <- load_immunity_coverages_by_age()
  coverage_infection_only <- coverage$infection * (1 - coverage$vaccination)
  coverage_vaccination_only <- coverage$vaccination * (1 - coverage$infection)
  coverage_infection_and_vaccination <- coverage$vaccination * coverage$infection

  infection_transmission_multiplier_by_age <- infection_ve_infection_only_by_age * coverage_infection_only +
    infection_ve_vaccination_only_by_age * coverage_vaccination_only +
    infection_ve_infection_and_vaccination_by_age * coverage_infection_and_vaccination

  onward_transmission_multiplier_by_age <- onward_ve_infection_only_by_age * coverage_infection_only +
    onward_ve_vaccination_only_by_age * coverage_vaccination_only +
    onward_ve_infection_and_vaccination_by_age * coverage_infection_and_vaccination

  # combine into a matrix to scale the NGM
  transmission_reduction_matrix <- infection_transmission_multiplier_by_age %*% t(onward_transmission_multiplier_by_age)

  # 3. use these average VEs by age with the NGM to get the expected reduction in
  # transmission due to immunity.

  # copy this functionality over from TP/Reff or conmat
  ngm_without_immunity <- load_ngm_without_immunity()

  # scale the NGM with average VEs for each age group
  ngm_with_immunity <- ngm_without_immunity * transmission_reduction_matrix

  # we could calibrate the NGM so that the reproduction number without the immunity added
  # is 1, so we don't need to recalculate this all the time, or just put the fixed value in here
  reproduction_number_without_immunity <- get_reproduction_number(ngm_without_immunity)

  # Compute reduction in dominant eigenvalue of the NGM
  reproduction_number_with_immunity <- get_reproduction_number(ngm_with_immunity)

  expected_reduction <- reproduction_number_with_immunity / reproduction_number_without_immunity

  # 4. define a likelihood between the expected reduction in TP due to immunity
  # and the distributions of these ratios from the TP/Reff model (posterior of
  # the ratios of the TP without immunity and corresponding reff estimate) in
  # the same time and state.

  # fill these in from posterior from TP/Reff model:
  # ratio_estimate <- 1
  # ratio_estimate <- 2
  distribution(ratio_estimate) <- normal(expected_reduction, ratio_estimate_sd)

  NULL

}
