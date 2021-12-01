#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param c50_vec
#' @param mean_log10_neut_vec
#' @param sd_log10_neut
#' @param log_k
#' @param lower
#' @param upper
#' @return
#' @author Nick Golding
#' @export
gaussian_ve_integrator <- function(
  c50_vec,
  mean_log10_neut_vec,
  sd_log10_neut,
  log_k,
  lower,
  upper
) {

  # dimensions and quadrature rules
  n_obs <- length(mean_log10_neut_vec)
  quads <- get_quad_rules(n_obs, lower = lower, upper = upper)
  n_quads <- length(quads$values)

  # expand out the vector parameters of the logit-normal density to matrices
  if (is.vector(c50_vec) & is.vector(mean_log10_neut_vec)) {
    c50_vec <- as.matrix(c50_vec)
    mean_log10_neut_vec <- as.matrix(mean_log10_neut_vec)
    vector_input <- TRUE
  } else {
    vector_input <- FALSE
  }

  repeater <- rep(1, n_quads)
  c50_mat <- c50_vec[, repeater]
  mean_log10_neut_mat <- mean_log10_neut_vec[, repeater]

  # and expand out the integration points to match
  values_matrix <- t(replicate(n_obs, quads$values))

  # get function values in matrix
  function_values <- logit_normal_density(
    x = values_matrix,
    c50 = c50_mat,
    mean_log10_neut = mean_log10_neut_mat,
    sd_log10_neut = sd_log10_neut,
    log_k = log_k
  )

  weights <- quads$weights

  # if we're doing this with greta arrays, we need to work around an issue with
  # greta checking matrix multiplly dimensions too early
  if (inherits(function_values, "greta_array")) {
    weights <- as_data(quads$weights)
  }

  ves <- function_values %*% weights

  if(vector_input) {
    dim(ves) <- NULL
  }

  ves

}
