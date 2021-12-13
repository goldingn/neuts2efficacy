#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ve_estimates
#' @return
#' @author Nick Golding
#' @export
# convert estimates into logits and compute logit sd
prep_ve_data_for_modelling <- function(ve_estimates) {

  ve_estimates %>%
    filter(
      variant == "delta",
      product != "mRNA booster"
    ) %>%
    mutate(
      # ensure VE upper bounds are not be 100% (rounding error)
      across(
        starts_with("ve"),
        ~pmin(.x, 0.99)
      ),
      # convert to logit scale
      across(
        starts_with("ve"),
        .fns = list(logit = qlogis)
      ),
      # compute the observation sd from 95% CIs
      ve_logit_sd = ((ve_upper_logit - ve_logit) + (ve_logit - ve_lower_logit)) / (2 * qnorm(0.975)),
      .after = ve_logit
    )

}
