#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param khoury_natmed_estimates
#' @return
#' @author Nick Golding
#' @export
get_neut_ratios_vaccine <- function() {

  cromer_url <- "https://github.com/InfectionAnalytics/SARS-CoV-2-Variants-and-Boosting---Lancet-Microbe/raw/main/Figure2-4andS4/SummaryTable_Efficacy_NeutRatio_SD_SEM.csv"
  cromer_file <- "data/SummaryTable_Efficacy_NeutRatio_SD_SEM_cromer.csv"

  if (!file.exists(cromer_file)) {
    download.file(cromer_url, cromer_file)
  }

  cromer_file %>%
    read_csv(
      col_types = cols(
        .default = col_skip(),
        Study = col_character(),
        SD = col_double(),
        NeutRatio_cens = col_double(),
        NumberIndividuals_Vaccine = col_double(),
        PooledSD = col_double()
      )
    ) %>%
    filter(
      Study %in% c("AstraZeneca", "Pfizer", "Convalescence")
    ) %>%
    mutate(
      immunity = case_when(
        Study == "AstraZeneca" ~ "AZ_dose_2",
        Study == "Pfizer" ~ "Pfizer_dose_2",
        Study == "Convalescence" ~ "infection"
      )
    ) %>%
    mutate(
      # get the standard error of the mean number of neuts *not of the ratio of convalescent*
      se = SD / sqrt(NumberIndividuals_Vaccine)
    ) %>%
    select(
      immunity,
      # mean of the normal distribution over log10 neuts, normalised to fold of
      # convalescent
      mean_log10_ratio_neut = NeutRatio_cens,
      # standard error of the mean of the normal distribution over log10 neuts
      sem_log10_neut = se,
      # shared SD of the normal distribution over log10 neuts
      sd_log10_ratio_neut = PooledSD
    ) %>%
    filter(
      immunity != "infection"
    )

}
