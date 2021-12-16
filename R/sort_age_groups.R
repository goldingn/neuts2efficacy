#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
sort_age_groups <- function(age_groups) {
  start <- str_split_fixed(age_groups, "-", 2)[, 1]
  order <- str_order(start, numeric = TRUE)
  age_groups[order]
}
