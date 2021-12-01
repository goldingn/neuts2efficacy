#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param x
#' @param mean
#' @param sd
#' @return
#' @author Nick Golding
#' @export
# greta function for the normal PDF
normal_density <- function(x, mean, sd) {

  op <- greta::.internals$nodes$constructors$op
  op(
    "normal_density",
    x, mean, sd,
    tf_operation = "tf_normal_density"
  )

}
