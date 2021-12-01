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

  khoury_natmed_estimates %>%
    select(
      product = Study,
      log10_ratio_neut = NeutRatio_cens
    ) %>%
    filter(
      product %in% c("AstraZeneca", "Pfizer", "Convalescence")
    ) %>%
    mutate(
      product = case_when(
        product == "AstraZeneca" ~ "AZ",
        product == "Pfizer" ~ "Pfizer",
        product == "Convalescence" ~ "Infection"
      ),
      ratio_neut = 10 ^ log10_ratio_neut
    ) %>%
    mutate(
      dose = if_else(product == "infection", 1, 2),
      .after = product
    )
}
