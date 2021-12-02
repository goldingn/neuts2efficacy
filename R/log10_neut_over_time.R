#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param time
#' @param maximum_log10_neut
#' @param decay
#' @return
#' @author Nick Golding
#' @export
# compute the log10 neut titre at a given time since the latest dose
log10_neut_over_time <- function (time, maximum_log10_neut, decay){
  # equivalent to: log10(10 ^ maximum_log10_neut * exp(-decay * time))
  maximum_log10_neut - decay * time / log(10)
}
