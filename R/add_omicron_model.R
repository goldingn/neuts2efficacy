#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param neut_model_initial
#' @return
#' @author Nick Golding
#' @export
# add Omicron model components to the basic model
add_omicron_model <- function(neut_model) {

  # get the omicron neut fold
  omicron_log10_neut_fold <- neut_model$model_objects$omicron_log10_neut_fold

  # define parameter for ratio of R0s (omicron / delta) with mode of 1
  R0_ratio <- normal(1, 1, truncation = c(0, Inf))

  # parameter for fraction with some immunity (from prior infection and/or vaccination)
  fraction_immune <- normal(0.8, 0.05, truncation = c(0, 1))

  # parameters for the increase in neut fold of average immune people in ZA, to
  # account for the potential boosting effect of large waves of different
  # variants (WT, beta, delta) and vaccination
  za_baseline_immunity_log10_neut_fold <- normal(0, 1, truncation = c(0, Inf))

  # add a parameter for the degree to which transmission in South Africa is
  # reduced by things other than immunity - e.g. hygiene and avoidance measures,
  # and isolation of cases
  za_distancing_effect <- normal(0.2, 0.1, truncation = c(0, 0.5))

  # correction factor for biases in the reinfection log hazard ratio -
  # informative prior to get the Delta wave reinfection hazard ratio in the same
  # ballpark as waned vaccine against symptomatic disease (AZ - Pfizer VEs 50 & 75% from
  # Andrews, hazard ratio 0.25-0.5)
  reinfection_correction_bounds <- log(0.1 / c(0.25, 0.5))
  reinfection_correction_estimate <- mean(reinfection_correction_bounds)
  reinfection_correction_sd <- mean(abs(reinfection_correction_bounds - reinfection_correction_estimate)) / 1.96
  reinfection_correction <- normal(reinfection_correction_estimate, reinfection_correction_sd)

  # define likelihood for R0 ratio based on ratio of secondary attack rates among
  # unvaccinated (assuming R0 is linear in the transmission probability)

  # Lyngse et al. https://www.medrxiv.org/content/10.1101/2021.12.27.21268278v1
  # give an estimated 'HSAR' ratio of 1.17 (0.99-1.38), but the way they
  # calculated this is equivalent to a HFAR, rather than a transmission
  # probability. So apply an approximate correction in the likelihood to convert
  # from a transmission probability ratio to an HFAR ratio. Base this on Fig 1
  # from Sharker et al. https://doi.org/10.1371/journal.pcbi.1008601, assuming
  # linearity in the SAR-FAR ratio over the range of FAR values.

  # get SARs for both variants
  sar_delta <- normal(0.3, 1, truncation = c(0, 1))
  sar_omicron <- sar_delta * R0_ratio

  # convert to FAR ratio (what Lyngse measured)
  far_delta <- far_to_sar(sar_delta)
  far_omicron <- far_to_sar(sar_omicron)
  far_ratio <- far_omicron / far_delta

  # define the likelihood over this data
  observed_far_ratio <- 1.17
  far_ratio_bounds <- c(0.99, 1.38)
  far_ratio_sd <- mean(abs(far_ratio_bounds - observed_far_ratio)) / 1.96

  distribution(observed_far_ratio) <- normal(far_ratio, far_ratio_sd, truncation = c(0, Inf))

  # # uncertainty in the log hazard ratio data
  # reinfection_log_hazard_ratio_sd <- normal(0, 1, truncation = c(0, Inf))

  # assume that immunity is largely driven by previous infection, since
  # vaccination rates are primarily in the older population, and vaccine uptake
  # is likely to be in the higher socioeconomic status individuals. This is
  # unlikely to affect the results however, since vaccine-derived immunity and
  # infection derived immunity are simmilar, and we are only looking at the
  # waning ratio anyway

  # define delay since peak for Omicron period, but integrate over a
  # distribution: mean of 120 days since peak (time since peak of previous
  # wave), but 95% within 80-160 days
  za_waning_days_raw <- normal(0, 1)
  za_waning_days <- 4 * 30 + za_waning_days_raw * 20
  # quantile(calculate(za_waning_days, nsim = 1e5)[[1]], c(0.025, 0.975))

  # get log10 neuts of convalescent (0) against omicron
  # (omicron_log10_neut_fold), this number of days post-peak
  omicron_log10_neuts <- log10_neut_over_time(
    time = za_waning_days,
    maximum_log10_neut = za_baseline_immunity_log10_neut_fold + omicron_log10_neut_fold,
    decay = neut_model$model_objects$neut_decay
  )

  delta_log10_neuts <- log10_neut_over_time(
    time = za_waning_days,
    maximum_log10_neut = za_baseline_immunity_log10_neut_fold,
    decay = neut_model$model_objects$neut_decay
  )

  # predict VEs for, for convalescent against omicron and against Delta for the
  # stated delays
  omicron_ves <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = rep(omicron_log10_neuts, 5),
    sd_log10_neut = neut_model$model_objects$sd_log10_neut_titres,
    log_k = neut_model$model_objects$log_k,
    c50_vec = neut_model$model_objects$c50s,
    method =  "gaussian"
  )

  delta_ves <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = rep(delta_log10_neuts, 5),
    sd_log10_neut = neut_model$model_objects$sd_log10_neut_titres,
    log_k = neut_model$model_objects$log_k,
    c50_vec = neut_model$model_objects$c50s,
    method =  "gaussian"
  )

  symptoms_idx <- match("symptoms", neut_model$lookups$outcome)
  acquisition_idx <- match("acquisition", neut_model$lookups$outcome)
  transmission_idx <- match("transmission", neut_model$lookups$outcome)

  # compute overall vaccine reduction in transmission from VEs for acquisition and transmisson
  # one minus the overall VE against acquisitoin and transmission
  omicron_vaccine_reduction <- prod(1 - omicron_ves[c(acquisition_idx, transmission_idx)])
  delta_vaccine_reduction <- prod(1 - delta_ves[c(acquisition_idx, transmission_idx)])

  # define likelihood on ratio of Reffs between Omicron and Delta in same
  # period, with expectation of:

  # immune_multiplier = fraction_immune * (1 - ve) + (1 - fraction_immune)
  # immune_multiplier_omicron = fraction_immune * (1 - ve_omicron) + (1 - fraction_immune)
  # immune_multiplier_delta = fraction_immune * (1 - ve_delta) + (1 - fraction_immune)

  # immune_multiplier_ratio = immune_multiplier_omicron / immune_multiplier_delta

  # Reff_omicron  =  R0_omicron * other_reductions * immune_multiplier_omicron
  # Reff_delta  =  R0_delta * other_reductions * immune_multiplier_delta

  # other reductions and fraction with immunity apply equally to both (assuming
  # GI  etc. the same), so:

  # Reff_ratio = Reff_omicron / Reff_delta
  #            = (R0_omicron * other_reductions * immune_multiplier_omicron) / (R0_delta * other_reductions * immune_multiplier_delta)
  #            = (R0_omicron * immune_multiplier_omicron) / (R0_delta * immune_multiplier_delta)
  #            = (R0_omicron / R0_delta) * (immune_multiplier_omicron / immune_multiplier_delta)
  #            = R0_ratio * immune_multiplier_ratio

  # where:
  # R0_ratio  = R0_omicron / R0_delta
  # ve_multiplier_ratio = (1 - ve_omicron) / (1 - ve_delta)

  # compute ascertainment reduction and transmission reduction separately and
  # combine, to mimic NGM approach
  fraction_nonimmune <- 1 - fraction_immune

  transmission_multiplier_omicron <- fraction_immune * (1 - omicron_ves[c(transmission_idx)]) + fraction_nonimmune
  acquisition_multiplier_omicron <- fraction_immune * (1 - omicron_ves[c(acquisition_idx)]) + fraction_nonimmune

  transmission_multiplier_delta <- fraction_immune * (1 - delta_ves[c(transmission_idx)]) + fraction_nonimmune
  acquisition_multiplier_delta <- fraction_immune * (1 - delta_ves[c(acquisition_idx)]) + fraction_nonimmune


  # get multiplier on probability of transmission per contact between an
  # infected person and an uninfected person (identical to ratio of NGM
  # eigenvalues in homogenously mixing population)
  immune_multiplier_omicron = transmission_multiplier_omicron * acquisition_multiplier_omicron
  immune_multiplier_delta = transmission_multiplier_delta * acquisition_multiplier_delta

  immune_multiplier_ratio <- immune_multiplier_omicron / immune_multiplier_delta

  expected_reff_ratio <- R0_ratio * immune_multiplier_ratio

  # define likelihood on ratio of Reff between Omicron and Delta

  # ratio of reffs
  # https://twitter.com/cap1024/status/1466840869852651529
  # https://drive.google.com/file/d/1hA6Mec2Gq3LGqTEOj35RqSeAb_SmXpbI/view
  # match most recent 95% CI
  reff_ratio_bounds <- c(1.3, 2.9)
  reff_ratio_estimate <- mean(reff_ratio_bounds)
  reff_ratio_sd <- mean(abs(reff_ratio_bounds - reff_ratio_estimate)) / 1.96
  # distribution(reff_ratio_estimate) <- normal(expected_reff_ratio, reff_ratio_sd)

  # define likelihood on reinfection hazard ratios with expectation of VE for
  # symptomatic disease (not acquisition, because of detection method), and with
  # a correction factor for biases in this estimate (reinfection hazard ratio
  # for Delta seems very low, compared with VEs).
  reinfection_log_hazard_ratio_omicron <- log(0.3)
  reinfection_log_hazard_ratio_delta <- log(0.1)
  # assume these estimates have low uncertainty, conditional on the model
  # parameters (we have a ratio parameter for uncertainty due to
  # misspecification)
  reinfection_log_hazard_ratio_sd <- 0.25

  omicron_expected_log_hazard_ratio <- log(1 - omicron_ves[symptoms_idx]) + reinfection_correction
  delta_expected_log_hazard_ratio <- log(1 - delta_ves[symptoms_idx]) + reinfection_correction

  # distribution(reinfection_log_hazard_ratio_omicron) <- normal(omicron_expected_log_hazard_ratio, reinfection_log_hazard_ratio_sd)
  # distribution(reinfection_log_hazard_ratio_delta) <- normal(delta_expected_log_hazard_ratio, reinfection_log_hazard_ratio_sd)

  # define a likelihood on Delta R0 vs Reff. Assume a Delta R0 of 6, with 20%
  # reduction (distancing, mask wearing) for ZA,
  R0_delta <- 6
  expected_log_reff_delta <- log(R0_delta * (1 - za_distancing_effect) * immune_multiplier_delta)
  log_reff_delta <- log(0.8)
  log_reff_delta_sd <- 0.2
  # distribution(log_reff_delta) <- normal(expected_log_reff_delta, log_reff_delta_sd)

  # redefine the greta model, tracing all the things that were traced in the previous version
  new_traced <- module(
    omicron_log10_neut_fold,
    R0_ratio,
    za_baseline_immunity_log10_neut_fold,
    za_waning_days,
    za_distancing_effect,
    fraction_immune,
    reinfection_correction,
    omicron_vaccine_reduction,
    delta_vaccine_reduction,
    expected_reff_ratio,
    omicron_expected_log_hazard_ratio,
    delta_expected_log_hazard_ratio,
    immune_multiplier_omicron,
    immune_multiplier_delta
  )

  old_traced_names <- names(neut_model$greta_model$dag$target_nodes)
  old_traced <- neut_model$model_objects[old_traced_names]
  traced <- c(old_traced, new_traced)

  # this is ugly, but the language-based version was too much effort
  call <- sprintf("greta::model(%s)", paste(names(traced), collapse = ", "))
  neut_model$greta_model <- eval(parse(text = call), envir = traced)

  neut_model$model_objects <- c(neut_model$model_objects, new_traced)

  neut_model
}
