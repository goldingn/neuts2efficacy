#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param khoury_natmed_estimates
#' @return
#' @author Nick Golding
#' @export
get_khoury_neut_ratios <- function(khoury_natmed_estimates) {

  khoury_natmed_file <- "data/SummaryTable_Efficacy_NeutRatio_SD_SEM.RDS"

  khoury_natmed_file %>%
    readRDS() %>%
    filter(
      Study %in% c("AstraZeneca", "Pfizer", "Convalescence")
    ) %>%
    mutate(
      Study = case_when(
        Study == "AstraZeneca" ~ "AZ",
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
    )

}
