test_that("log10_neut_over_time works", {

  n <- 10
  time <- seq(0, 100, length.out = n)
  max_log10_neut <- rnorm(n)
  halflife <- rnorm(n, 108, 10)
  decay <- log(2) / halflife

  # calculate current neuts the long way
  max_neut <- 10 ^ max_log10_neut
  current_neut <- max_neut * exp(-decay * time)
  current_log10_neut_long <- log10(current_neut)

  # and the short and computationally efficient way in the function
  current_log10_neut_short <- log10_neut_over_time(
    time = time,
    maximum_log10_neut = max_log10_neut,
    decay = decay
  )

  expect_equal(current_log10_neut_short, current_log10_neut_long)


})
