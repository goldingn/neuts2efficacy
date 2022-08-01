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
# an analytic and differentiable version of log10_neut_over_time_varying_numeric
# - both implement a waning function on the log10 neuts were therate of decay
# plateaus
log10_neut_over_time_varying_analytic <- function (time,
                                                   maximum_log10_neut,
                                                   max_decay,
                                                   min_decay,
                                                   end_decline = 365,
                                                   start_decline = 250){

  # compute the analytic solution to the decay rate in the initial, declining,
  # and final stages in a differentiable and vectorised way, to make greta happy

  # the analytic solution to exponential decay:
  #   dN/dt = -lambda * N
  # is:
  #   N[t] = N[0] * exp(-lambda * t))
  # The analytic solution to the same function but with a linearly varying decay
  # rate:
  #   lambda[t] = a + b * t
  #   dN/dt = -lambda[t] * N
  #   dN/dt = -(a + b * t) * N
  # is:
  #   N[t] = N[0] * exp(-a * t + -0.5 * b * t^2)
  # Tas Symons solved this, following the same method as described here:
  # https://en.wikipedia.org/wiki/Exponential_decay#Solution_of_the_differential_equation

  # the following code gets the neuts for the initial constant decay phase (part
  # 1), then in the declining decay phase (part 2) and then for the final phase
  # (part 3); these are then spliced together using differentiable masking
  # variables. This is a pain in the arse, but it's done because: a) the integer
  # step method in log10_neut_over_time_varying_numeric is biased down in the
  # final phase, for the Khoury et al.'s Figure 3b parameters (only by about
  # 0.001, but I'm fussy about these things) b) ODE solvers are horribly slow,
  # and c) we can't use Boolean logic in the likelihood of HMC (it causes a
  # non-differentiable likelihood), and we want to be able
  # to learn the decay plateau parameters.

  # for the first constant decay part:
  # equivalent to: log10(10 ^ maximum_log10_neut * exp(-decay * time))
  part1_time <- time
  log10_neut_part1 <- maximum_log10_neut - max_decay * part1_time / log(10)

  # for the second, linearly declining decay part, starting at the end of the next phase
  log10_neut_part2_start <- maximum_log10_neut - max_decay * start_decline / log(10)

  # get parameters
  decline_duration <- end_decline - start_decline
  linear_decay_slope <- (min_decay - max_decay) / decline_duration
  linear_decay_intercept <- max_decay

  # get neuts during this phase
  part2_time <- time - start_decline
  # equivalent to:
  # neut_part2_start <- 10^log10_neut_part2_start
  # neut_part2 <- neut_part2_start * exp(-linear_decay_intercept * part2_time + -0.5 * linear_decay_slope * part2_time ^ 2)
  # log10_neut_part2 <- log10(neut_part2)
  log10_neut_part2 <- log10_neut_part2_start + (-linear_decay_intercept * part2_time + -0.5 * linear_decay_slope * part2_time ^ 2) / log(10)

  # for the final, constant phase
  # equivalent to:
  # neut_part3_start <- neut_part2_start * exp(-linear_decay_intercept * decline_duration + -0.5 * linear_decay_slope * decline_duration ^ 2)
  # log10_neut_part3_start <- log10(neut_part3_start)
  log10_neut_part3_start <- log10_neut_part2_start + (-linear_decay_intercept * decline_duration + -0.5 * linear_decay_slope * decline_duration ^ 2) / log(10)

  part3_time <- time - end_decline
  log10_neut_part3 <- log10_neut_part3_start - min_decay * part3_time / log(10)

  # now get the correct part, with differentiable masks (ifelse etc. are not
  # differentiable, so cannot be used in HMC where end_decline is an unknown
  # parameters)
  mask_1 <- ilogit((start_decline - time) * 1e3)
  mask_3 <- ilogit((time - end_decline) * 1e3)
  mask_2 <- 1 - (mask_1 + mask_3)

  log10_neut <- log10_neut_part1 * mask_1 +
    log10_neut_part2 * mask_2 +
    log10_neut_part3 * mask_3

  log10_neut

}
