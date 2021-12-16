#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param nameme1
#' @return
#' @author Nick Golding
#' @export
get_ngm <- function() {

  if (file.exists("outputs/australia_ngm.RDS")) {
    australia_ngm <- readRDS("outputs/australia_ngm.RDS")
  } else {
    australia_ngm <- make_ngm()
    saveRDS(australia_ngm, "outputs/australia_ngm.RDS")
  }

  australia_ngm

}
