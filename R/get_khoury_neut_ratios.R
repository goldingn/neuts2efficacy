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
      mean_log10_ratio_neut = NeutRatio_cens,
      se_mean_log10_ratio_neut = SEM,
      sd_log10_ratio_neut = PooledSD
    ) %>%
    filter(
      product %in% c("AstraZeneca", "Pfizer")
    ) %>%
    mutate(
      product = case_when(
        product == "AstraZeneca" ~ "AZ",
        product == "Pfizer" ~ "Pfizer"
      )
    ) %>%
    mutate(
      dose = if_else(product == "infection", 1, 2),
      .after = product
    )
}
