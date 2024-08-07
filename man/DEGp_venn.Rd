% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deg.plot.R
\name{DEGp_venn}
\alias{DEGp_venn}
\title{Venn Plot}
\usage{
DEGp_venn(datalist, transparency = 0.6)
}
\arguments{
\item{datalist}{list, such as list(diffan$deg$symbol, diffagek, ...), no more than 5-dimensional-list}

\item{transparency}{number, transparency, 0-1, also can be NULL, transparency = NULL means no fill color}
}
\value{
plot
}
\description{
Plot Venn 2-5 group data
}
\examples{
\dontrun{
lzpv <- list(a=1:10, b=5:16, c=15:23, d=c(1,2,4,8,8,16:19,30:34))
p <- lzplot_veen(datalist = lzpv, transparency = 0.6)
}
}
\author{
Jiang
}
