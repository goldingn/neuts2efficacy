test_that("integration methods agree on R objects", {

  n <- 10
  c50s <- rnorm(n)
  neuts <- rnorm(10, c50s, 0.5)

  adaptive <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = neuts,
    c50_vec = c50s,
    sd_log10_neut = 0.5,
    log_k = 0,
    method = "adaptive"
  )

  gaussian <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = neuts,
    c50_vec = c50s,
    sd_log10_neut = 0.5,
    log_k = 0,
    method = "gaussian"
  )

  expect_equal(gaussian, adaptive)

})

test_that("gaussian integration agrees between R and greta", {

  n <- 10
  c50s <- rnorm(n)
  neuts <- rnorm(10, c50s, 0.5)

  r_result <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = neuts,
    c50_vec = c50s,
    sd_log10_neut = 0.5,
    log_k = 0,
    method = "adaptive"
  )

  greta_result_ga <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = as_data(neuts),
    c50_vec = as_data(c50s),
    sd_log10_neut = as_data(0.5),
    log_k = as_data(0),
    method = "gaussian"
  )

  greta_result <- calculate(greta_result_ga)[[1]]
  dim(greta_result) <- NULL

  expect_equal(greta_result, r_result, tolerance = 1e-6)

})
