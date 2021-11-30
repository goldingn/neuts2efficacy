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
log10_neut_over_time = function (time, maximum_log10_neut, decay){
  maximum_neut <- 10 ^ maximum_log10_neut
  current_neut <- maximum_neut * exp(-decay * time)
  current_log10_neut <- log10(current_neut)
  current_log10_neut
}
