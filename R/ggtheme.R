#' Update ggplot theme text
#' @description Update ggplot theme text
#'
#' @param size number font size, default is 14
#' @param color character font color, default is "black", can be rgb value, such as "#00FFFF"
#' @param family character font family, default is "serif"
#' @importFrom ggplot2 theme
#'
#' @return ggtheme
#' @export
#' @author Jiang
ggtheme.update.text <- function(size = 14, color = "black", family = "serif") {
  list(
    theme(
      axis.text = element_text(size = size, color = color, family = family),
      title = element_text(size = size, face = "bold", family = family))
  )
}

#' Update ggplot theme text legend
#' @description Update ggplot theme text
#'
#' @param size number font size, default is 8
#' @param color character font color, default is "black", can be rgb value, such as "#00FFFF"
#' @param family character font family, default is "serif"
#' @importFrom ggplot2 theme
#'
#' @return ggtheme
#' @export
#' @author Jiang
ggtheme.update.text.legend <- function(size = 8, color = "black", family = "serif") {
  list(
    theme(
      legend.text = element_text(size = size, color = color, family = family)
    )
  )
}