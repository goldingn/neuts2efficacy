#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
get_ve_estimates <- function() {

  bind_rows(
    get_ve_estimates_andrews_delta(),
    get_ve_estimates_pouwels_delta(),
    get_ve_estimates_eyre_delta(),
    get_ve_estimates_andrews_omicron()
  ) %>%
    rowwise() %>%
    mutate(
      days = mean(days_earliest:days_latest),
      .after = dose
    ) %>%
    ungroup() %>%
    filter(
      # keep only those at least 14 days post vaccination
      days >= 0,
      # remove those with negative point estimates
      ve >= 0,
      # remove dose 2 estimates from the Andrews Omicron paper, since these will
      # duplicate the earlier Delta paper
      !(source == "andrews_omicron" & variant == "delta" & dose < 3),
      # remove the AZ estimates from the Andrews Omicron paper, since these will
      # be biased towards the more vulnerable (older, more co-morbidities)
      # population and have very small sample sizes
      !(source == "andrews_omicron" & product == "AZ")
    ) %>%
    mutate(
      # make all VE lower bounds positive at lowest reposrted resolution
      ve_lower = pmax(ve_lower, 0.001),
      # recode boosters as the same product (assume same VEs)
      product = case_when(
        dose == 3 ~ "mRNA booster",
        TRUE ~ product
      )
    )

}
