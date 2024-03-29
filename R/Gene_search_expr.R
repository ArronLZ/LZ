#' Plot Gene Expression Density
#' @description Plot Gene Expression Density
#' 
#' @param df data.frame, the data.frame of the clinical and gene expression data
#' @param tagetGene character, the gene name in the data.frame
#' @import ggplot2
#' @include ggtheme.R
#'
#' @return ggplot2 object
#' @export
#'
#' @examples #
Gene_dens <- function(df, tagetGene) {
  five <- fivenum(df[, tagetGene])
  ggplot(df, aes(x = df[, tagetGene])) +
    geom_density(fill = "gray", alpha = 0.5) +
    geom_vline(xintercept = five, col = "blue", size=1) +
    geom_vline(xintercept = median(df[, tagetGene]), col = "red", size=1) +
    geom_vline(xintercept = mean(df[, tagetGene]), col = "yellow", size=1, lty=2) +
    annotate("text", x = five, y = 0, 
             label = c(paste0("Min:", round(five[1],2)), 
                       paste0("Q1:", round(five[2],2)), 
                       paste0("Median:", round(five[3],2)), 
                       paste0("Q3:", round(five[4],2)), 
                       paste0("Max:", round(five[5],2)) ), 
             color = c("black", "black", "red", "black", "black"),
             hjust = 1, vjust = c(-1,-2.5,-4.5,-6,-8)) +
    annotate("text", x = mean(df[, tagetGene]), y = 0, 
             label = paste0("Mean:", round(mean(df[, tagetGene]))), color="yellow",
             hjust = 1, vjust = -1) +
    mytheme_prism() + xlab(label = tagetGene) + ylab(label = NULL)
}






