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

  # define parameter for fold of neuts relative to main model
  omicron_log10_neut_fold <- normal(0, 10)

  # define parameter for ratio of R0s: mode of 1
  R0_ratio <- normal(1, 1, truncation = c(0, Inf))

  # correction factor for biases in the reinfection hazard ratio
  reinfection_correction <- normal(1, 1, truncation = c(0, Inf))

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
    maximum_log10_neut = 0 + omicron_log10_neut_fold,
    decay = neut_model$model_objects$neut_decay
  )

  delta_log10_neuts <- log10_neut_over_time(
    time = za_waning_days,
    maximum_log10_neut = 0,
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

  # Reff_delta  =  R0_delta * (1 - ve_delta) * fraction_with_immunity * other_reductions
  # Reff_omicron  =  R0_omicron * (1 - ve_omicron) * fraction_with_immunity * other_reductions

  # other reductions and fraction with immunity apply equally to both (assuming
  # GI  etc. the same), so:

  # Reff_ratio = Reff_omicron / Reff_delta
  #            = (R0_omicron * (1 - ve_omicron)) / (R0_delta * (1 - ve_delta))
  #            = (R0_omicron / R0_delta) * ((1 - ve_omicron) / (1 - ve_delta))
  #            = R0_ratio * ve_multiplier_ratio

  # where:
  # R0_ratio  = R0_omicron / R0_delta
  # ve_multiplier_ratio = (1 - ve_omicron) / (1 - ve_delta)

  ve_multiplier_ratio <- omicron_vaccine_reduction / delta_vaccine_reduction
  expected_reff_ratio <- R0_ratio * ve_multiplier_ratio

  # define likelihood on ratio of Reff between Omicron and Delta

  # ratio of reffs
  # https://twitter.com/cap1024/status/1466840869852651529
  # https://drive.google.com/file/d/1hA6Mec2Gq3LGqTEOj35RqSeAb_SmXpbI/view
  # match most recent 95% CI
  reff_ratio_bounds <- c(1.4, 3.2)
  reff_ratio_estimate <- mean(reff_ratio_bounds)
  reff_ratio_sd <- mean(abs(reff_ratio_bounds - reff_ratio_estimate)) / 1.96

  distribution(reff_ratio_estimate) <- normal(expected_reff_ratio, reff_ratio_sd)

  # define likelihood on reinfection hazard ratios with expectation of VE for
  # symptomatic disease (not acquisition, because of detection method), and with
  # a correction factor for biases in this estimate (reinfection hazard ratio
  # for Delta seems very low, compared with VEs).
  reinfection_hazard_ratio_omicron <- 0.3
  reinfection_hazard_ratio_delta <- 0.1
  # set the observation SD such that there is a 1 in 1000 chance that the reinfection difference is completely spurious and there's no real difference
  reinfection_hazard_ratio_sd <- (reinfection_hazard_ratio_omicron - reinfection_hazard_ratio_delta) / qnorm(1 - 1 / 1000)
  # 1 / (1 - pnorm(reinfection_hazard_ratio_omicron, reinfection_hazard_ratio_delta, reinfection_hazard_ratio_sd))
  # hist(rnorm(1e5, reinfection_hazard_ratio_omicron, reinfection_hazard_ratio_sd), breaks = 100)
  # abline(v = c(reinfection_hazard_ratio_omicron, reinfection_hazard_ratio_delta))

  omicron_expected_hazard_ratio <- omicron_ves[symptoms_idx] * reinfection_correction
  delta_expected_hazard_ratio <- delta_ves[symptoms_idx] * reinfection_correction

  distribution(reinfection_hazard_ratio_omicron) <- normal(omicron_expected_hazard_ratio, reinfection_hazard_ratio_sd)
  distribution(reinfection_hazard_ratio_delta) <- normal(delta_expected_hazard_ratio, reinfection_hazard_ratio_sd)

  new_traced <- module(
    omicron_log10_neut_fold,
    R0_ratio,
    za_waning_days,
    reinfection_correction,
    expected_reff_ratio,
    omicron_expected_hazard_ratio,
    delta_expected_hazard_ratio
  )

  old_traced_names <- names(neut_model$greta_model$dag$target_nodes)
  old_traced <- neut_model$model_objects[old_traced_names]

  traced <- c(old_traced, new_traced)
  traced_names <- names(traced)
  call <- sprintf("greta::model(%s)",
                  paste(names(traced), collapse = ", "))

  neut_model$greta_model <- eval(parse(text = call), envir = traced)

  neut_model$model_objects <- c(neut_model$model_objects, new_traced)

  neut_model
}
