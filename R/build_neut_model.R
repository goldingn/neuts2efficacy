#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
# build a greta model to infer mean neuts, c50s, and the slope of the
# relationship between VEs and neut titres
build_neut_model <- function(ve_estimates, neut_ratios_vaccine) {

  # set up ve data for modelling
  ve_data_modelling <- ve_estimates %>%
    prep_ve_data_for_modelling()

  # record unique values of levels; parameter orders correspond to these
  lookups <- list(
    variant = unique(ve_data_modelling$variant),
    outcome = unique(ve_data_modelling$outcome),
    immunity = unique(ve_data_modelling$immunity)
  )

  # get indices to various objects to get everything in the correct order for
  # prediction
  indices <- list(
    variant_idx = match(ve_data_modelling$variant, lookups$variant),
    outcomes_idx = match(ve_data_modelling$outcome, lookups$outcome),
    immunity_idx = match(ve_data_modelling$immunity, lookups$immunity),
    neut_ratio_vaccine_idx = match(neut_ratios_vaccine$immunity, lookups$immunity)
  )

  # define greta model

  # halflife of the neutralising antibody titre. Khoury found 108 days from
  # convalescent this was initialising a bit strangely when defined as N(108, 10),
  # so define parameter as standard normmal and transform
  neut_halflife_raw <- normal(0, 1)
  neut_halflife <- 108 + neut_halflife_raw * 10

  # convert to exponential rate parameter
  neut_decay <- log(2) / neut_halflife

  # log slope of the logit mapping from neuts to VEs
  log_k <- normal(0, 1)

  # a shared additional observation error on logit VEs (since uncertainty in
  # estimates doesn't reflect differences in outcome definitions between studies
  # etc.)
  ve_logit_obs_sd_shared <- normal(0, 1, truncation = c(0, Inf))

  # c50s for different outcomes
  c50s <- normal(0, 1, dim = 5)

  # a fold multiplier for mRNA booster relative to Pfizer dose 2
  # (mean 5-fold, 95% CI 3 to 7-fold)
  booster_multiplier <- normal(5, 1, truncation = c(0, Inf))

  # a fold multiplier for neuts against Omicron, relative to Dellta. Expect this
  # to be negative (less effect against Omicron), but don't constrain it.
  omicron_log10_neut_fold <- normal(0, 1)

  # combine into vector
  delta_log10_neut_fold <- zeros(1)

  variant_log10_neut_fold_unordered <- c(omicron_log10_neut_fold, delta_log10_neut_fold)
  variant_order <- match(c("omicron", "delta"), lookups$variant)
  variant_log10_neut_fold <- variant_log10_neut_fold_unordered[variant_order]

  # mean log10 peak neuts for second doses of AZ and Pfizer use the Khoury mean
  # (of the log10 ratio to convalescent) and standard error (of the ratio to
  # convalescent)
  dose_2_mean_log10_neuts <- normal(
    mean = neut_ratios_vaccine$mean_log10_ratio_neut[indices$neut_ratio_vaccine_idx],
    sd = neut_ratios_vaccine$sem_log10_neut[indices$neut_ratio_vaccine_idx]
  )

  # difference in peak neut titre between first and second doses of each
  # (constrain dose 2 to not be lower than dose 1)
  dose_2_vs_1_mean_log10_neut_increase <- normal(0, 1, dim = 2, truncation = c(0, Inf))

  # mean log10 peak neuts for first doses of AZ and Pfizer respectively
  dose_1_mean_log10_neuts <- dose_2_mean_log10_neuts - dose_2_vs_1_mean_log10_neut_increase

  # and for a booster
  booster_mean_log10_neuts <- dose_2_mean_log10_neuts[2] + log10(booster_multiplier)

  # define these doses as AZ first, then Pfizer, and ensure the order is correct when comparing to data
  peak_mean_log10_neuts_unordered <- c(dose_1_mean_log10_neuts, dose_2_mean_log10_neuts, booster_mean_log10_neuts)
  neuts_order <- match(c("AZ_dose_1", "Pfizer_dose_1", "AZ_dose_2", "Pfizer_dose_2", "mRNA_booster"), lookups$immunity)
  peak_mean_log10_neuts <- peak_mean_log10_neuts_unordered[neuts_order]

  # pull out vectors of parameters

  # expand out the c50 parameters to match the data
  c50_vec <- c50s[indices$outcomes_idx]

  # get peak mean log10 neut titres for each vaccine and variant
  peak_mean_log10_neut_vec <- peak_mean_log10_neuts[indices$immunity_idx] +
    variant_log10_neut_fold[indices$variant_idx]

  # compute expected values on given days post peak
  mean_log10_neut_vec <- log10_neut_over_time(
    time = ve_data_modelling$days,
    maximum_log10_neut = peak_mean_log10_neut_vec,
    decay = neut_decay
  )

  # population standard deviation of log10 neut titres
  sd_log10_neut_titres <- neut_ratios_vaccine$sd_log10_ratio_neut[1]

  # expected VEs for all combinations
  ve_expected <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = mean_log10_neut_vec,
    c50_vec = c50_vec,
    sd_log10_neut = sd_log10_neut_titres,
    log_k = log_k,
    method = "gaussian"
  )

  # convert to logit scale
  ve_expected_logit <- log(ve_expected / (1 - ve_expected))

  # get the observation standard deviation for each observation; combining error
  # from both the estimation (known for each estimate) with error between studies
  # in outcome definitions and sample biases etc
  ve_logit_obs_sd <- sqrt(ve_data_modelling$ve_logit_sd ^ 2 + ve_logit_obs_sd_shared ^ 2)

  # likelihood for VEs, with observation error given by estimate SD
  distribution(ve_data_modelling$ve_logit) <- normal(ve_expected_logit, ve_logit_obs_sd)

  # define greta model object
  greta_model <- greta::model(
    omicron_log10_neut_fold,
    booster_multiplier,
    log_k,
    neut_halflife,
    c50s,
    dose_1_mean_log10_neuts,
    dose_2_mean_log10_neuts
  )

  # group together the model objects
  model_objects <- module(
    omicron_log10_neut_fold,
    variant_log10_neut_fold,
    booster_multiplier,
    neut_halflife,
    neut_decay,
    log_k,
    ve_logit_obs_sd_shared,
    c50s,
    sd_log10_neut_titres,
    dose_2_mean_log10_neuts,
    dose_1_mean_log10_neuts,
    peak_mean_log10_neuts,
    ve_expected,
    ve_logit_obs_sd
  )

  # group together greta arrays
  neut_model <- module(
    ve_data_modelling,
    lookups,
    indices,
    greta_model,
    model_objects
  )

}
