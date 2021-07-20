#' The indicator function of an event
#'
#' @param bool event or logical statement to indicate
#'
#' @description {Returns zero if the event does not occur and 1 if it does.}
#' @return numeric
#' @export indicator
indicator <- function(bool)
{
  ifelse(bool, 1, 0)
}
