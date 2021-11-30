#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param .
#' @return
#' @author Nick Golding
#' @export
# given data on waning of different VEs form Andrews et al., with C50s and
# initial neuts attached as columns, estimate the rate of decay of neuts by
# least squares, accounting for potential biases in their VEs (set mean ratio of
# observed and predicted VEs to be 0 for each outcome - equivalent to
# multiplicative correction factor)
find_decay_rate <- function(data) {

  optimiser_fun <- function(log_decay_rate, data) {

    data %>%
      mutate(
        decay = exp(log_decay_rate),
        current_log10_neut = log10_neut_over_time(
          time = days,
          maximum_log10_neut = optimal_log10_neut,
          decay = decay
        )
      ) %>%
      rowwise() %>%
      mutate(
        predicted_ve = ve_from_mean_log10_neut(
          mean_log10_neut = current_log10_neut,
          sd_log10_neut = sd_log10_neut_titres,
          log_k = log_k,
          c50 = c50
        )
      ) %>%
      ungroup() %>%
      group_by(
        outcome
      ) %>%
      # compute ratio to correct VEs for each outcome type
      mutate(
        correction = mean(ve / predicted_ve),
        predicted_ve = predicted_ve * correction
      ) %>%
      ungroup() %>%
      mutate(
        error = ve - predicted_ve
      ) %>%
      summarise(
        sum_of_squares = sum(error ^ 2)
      ) %>%
      pull(
        sum_of_squares
      )

  }

  result <- optimise(optimiser_fun, interval = c(-6, 0), data = data)
  exp(result$minimum)

}
