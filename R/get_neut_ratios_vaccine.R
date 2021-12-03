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

  khoury_natmed_url <- "https://github.com/InfectionAnalytics/COVID19-ProtectiveThreshold/raw/main/Comparing%20Severe%20vs%20Mild%20Infection/SummaryTable_Efficacy_NeutRatio_SD_SEM.csv"
  khoury_natmed_file <- "data/SummaryTable_Efficacy_NeutRatio_SD_SEM.csv"

  if (!file.exists(khoury_natmed_file)) {
    download.file(khoury_natmed_url, khoury_natmed_file)
  }

  khoury_natmed_file %>%
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
      Study %in% c("Astra", "Pfizer", "Convalescence")
    ) %>%
    mutate(
      Study = case_when(
        Study == "Astra" ~ "AZ",
        Study == "Pfizer" ~ "Pfizer",
        Study == "Convalescence" ~ "infection"
      )
    ) %>%
    mutate(
      dose = if_else(Study == "infection", 1, 2),
      .after = Study
    ) %>%
    mutate(
      # get the standard error of the mean number of neuts *not of the ratio of convalescent*
      se = SD / sqrt(NumberIndividuals_Vaccine)
    ) %>%
    select(
      product = Study,
      # mean of the normal distribution over log10 neuts, normalised to fold of
      # convalescent
      mean_log10_ratio_neut = NeutRatio_cens,
      # standard error of the mean of the normal distribution over log10 neuts
      sem_log10_neut = se,
      # shared SD of the normal distribution over log10 neuts
      sd_log10_ratio_neut = PooledSD
    ) %>%
    filter(
      product != "infection"
    )

}
