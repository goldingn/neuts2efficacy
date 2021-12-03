# checking the censored SD estimate is unbiased
est_params_censored_normal <- function() {

  parameters <- c("mean", "sd")
  truth <- c(runif(1, -3, 3), runif(1, 0.5, 1.5))

  x <- rnorm(50, truth[1], truth[2])
  lod <- quantile(x, 0.1)
  low <- x < lod
  x[low] <- lod

  est_naive <- c(mean(x), sd(x))

  censored_llik <- function(theta) {
    sum(pnorm(x[low], theta[1], theta[2], log.p = TRUE))
  }

  uncensored_llik <- function(theta) {
    sum(dnorm(x[!low], theta[1], theta[2], log = TRUE))
  }

  nll <- function(theta) {
    -(uncensored_llik(theta) + censored_llik(theta))
  }

  o <- nlm(
    f = nll,
    p = est_naive
  )

  est_mle <- o$estimate

  tibble(
    parameters,
    truth,
    est_mle,
    est_naive
  )

}

sims <- replicate(
  1000,
  est_params_censored_normal(),
  simplify = FALSE
)

sims %>%
  do.call(rbind, .) %>%
  filter(
    parameters == "sd"
  ) %>%
  mutate(
    error_mle = est_mle - truth,
    error_naive = est_naive - truth
  ) %>%
  pivot_longer(
    cols = starts_with("error"),
    names_to = "type",
    names_prefix = "error_",
    values_to = "error"
  ) %>%
  ggplot(
    aes(
      x = type,
      y = error
    )
  ) +
  geom_jitter(
    colour = grey(0.8)
  ) +
  stat_summary(
    aes(ymax = ..y.., ymin = ..y..),
    geom = "errorbar",
    fun = "mean",
    width = 0.75,
    size = 1
  ) +
  theme_minimal()



