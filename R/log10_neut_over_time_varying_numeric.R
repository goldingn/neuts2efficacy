#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param time
#' @param maximum_log10_neut
#' @param max_decay
#' @param min_decay
#' @param end_decline
#' @param start_decline
#' @return
#' @author Nick Golding
#' @export
# compute the log10 neuts with a time-varying decay - using the Khoury et al.
# decay model
log10_neut_over_time_varying_numeric <- function (time,
                                                  maximum_log10_neut,
                                                  max_decay,
                                                  min_decay,
                                                  end_decline = 365,
                                                  start_decline = 250){
  # construct a timeseries of decay rates from decay parameters, and compute
  # cumulative decay
  times <- seq_len(max(time))
  # construct a timeseries of decay rates from the decay parameters
  decay_timeseries <- decay_over_time(time = times,
                                      max_decay = max_decay,
                                      min_decay = min_decay,
                                      end_decline = end_decline,
                                      start_decline = start_decline)
  # compute the cumulative amount of decay at all times, anc convert to the
  # log10 neut titres at all times
  log10_neut_timeseries <- maximum_log10_neut + cumsum(log10(1 - decay_timeseries))
  # this^ is a numerically more stable version of this:
  # log10_neut_timeseries <- log10(10 ^ maximum_log10_neut * cumprod(1 - decay_timeseries))
  log10_neut_timeseries[time]
}

# parameterise a time-varying rate of decay, according to Khoury et al. Nature
# Med 2021: "This model uses the estimated half-life of SARS-CoV-2
# neutralization titer in convalescent subjects of 108d over the first 250d,
# after which the decay decreases linearly until it achieves a 10-year half-life
# (consistent with the long-term stability of antibody responses seen after
# other vaccines)"
decay_over_time <- function(time, max_decay, min_decay, end_decline = 365, start_decline = 250) {

  # calculate te end of the decline phase
  decline_duration <- end_decline - start_decline

  # define a linear function so that it is max_decay at start_linear, and
  # min_decay at end_linear, then use pmin and pmax to clamp it

  # calculate slope of linear decline part of decay curve
  linear_decay_slope <- (min_decay - max_decay) / decline_duration

  # calculate intercept of linear decline part of decay curve, extrapolating
  # back from the start
  linear_decay_intercept <- max_decay - linear_decay_slope * start_decline

  # combine to get the linear component
  linear_decay <- linear_decay_intercept + linear_decay_slope * time

  # clamp the decay over time to be constant before and after the linear period
  decay <- pmin(max_decay, pmax(min_decay, linear_decay))

  decay

}
