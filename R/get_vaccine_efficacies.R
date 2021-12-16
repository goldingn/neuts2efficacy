#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param vaccine_cohorts_now
#' @return
#' @author Nick Golding
#' @export
get_vaccine_efficacies <- function(vaccine_cohorts) {

  # load omicron parameters in wide format and subset to different parameter sets
  params_wide <- get_omicron_params_wide()

  neut_params_wide <- params_wide %>%
    select(
      omicron_scenario,
      log10_mean_neut_AZ_dose_1,
      log10_mean_neut_AZ_dose_2,
      log10_mean_neut_Pfizer_dose_1,
      log10_mean_neut_Pfizer_dose_2,
      log10_mean_neut_mRNA_booster,
      neut_decay
    )

  ve_params_wide <- params_wide %>%
    select(
      omicron_scenario,
      starts_with("c50"),
      log_k,
      sd_log10_neut_titres,
      omicron_log10_neut_fold
    )

  # compute the average neutralisation level (mean log10 neut fold of WT
  # convalescent) in each age group, scenario, and omicron scenario
  mean_neuts <- vaccine_cohorts %>%
    filter(
      !is.na(immunity)
    ) %>%
    full_join(
      tibble(
        omicron_scenario = c(
          "intermediate",
          "optimistic",
          "pessimistic"
        )
      ),
      by = character()
    ) %>%
    left_join(
      neut_params_wide,
      by = "omicron_scenario"
    ) %>%
    # compute the peak and waned log10 mean neuts for each cohort
    mutate(
      peak_neuts = case_when(
        immunity == "AZ_dose_1" ~ log10_mean_neut_AZ_dose_1,
        immunity == "AZ_dose_2" ~ log10_mean_neut_AZ_dose_2,
        immunity == "Pf_dose_1" ~ log10_mean_neut_Pfizer_dose_1,
        immunity == "Pf_dose_2" ~ log10_mean_neut_Pfizer_dose_2,
        immunity == "mRNA_booster" ~ log10_mean_neut_mRNA_booster
      )
    ) %>%
    mutate(
      neuts = log10_neut_over_time(
        time = days_ago,
        maximum_log10_neut = peak_neuts,
        decay = neut_decay
      )
    ) %>%
    select(
      -starts_with("log10_mean_neut"),
      -peak_neuts,
      -neut_decay
    ) %>%
    # average the mean neuts over cohorts and scenarios
    group_by(
      scenario, omicron_scenario, age_band
    ) %>%
    summarise(
      neuts = weighted.mean(neuts, num_people),
      .groups = "drop"
    )

  # now compute VEs against each outcome, for Omicron and Delta
  ves <- mean_neuts %>%
    left_join(
      ve_params_wide,
      by = "omicron_scenario"
    ) %>%
    # for omicron, adjust down the neuts
    full_join(
      tibble(
        variant = c("Delta", "Omicron")
      ),
      by = character()
    ) %>%
    mutate(
      neuts = case_when(
        variant == "Omicron" ~ neuts + omicron_log10_neut_fold,
        TRUE ~ neuts
      )
    ) %>%
    # compute all the VEs in one shot with Gaussian integration
    pivot_longer(
      cols = starts_with("c50"),
      names_to = "outcome",
      values_to = "c50",
      names_prefix = "c50_"
    ) %>%
    mutate(
      ve = ve_from_mean_log10_neut(
        mean_log10_neut_vec = neuts,
        sd_log10_neut = sd_log10_neut_titres,
        log_k = log_k,
        c50_vec = c50,
        method = "gaussian"
      )
    ) %>%
    select(
      -neuts,
      -log_k,
      -sd_log10_neut_titres,
      -omicron_log10_neut_fold,
      -c50
    )

  ves

}
