% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surv_km.R
\name{Surv_diffG}
\alias{Surv_diffG}
\title{Survival analysis and plot}
\usage{
Surv_diffG(
  df,
  col_name,
  num.tran = 365,
  break.x = 2.5,
  main.text,
  color = c("#eaa8a8", "#729fc9")
)
}
\arguments{
\item{df}{data.frame, the data.frame of the clinical and gene expression data}

\item{col_name}{character, the column name of the group column, default is "Group"}

\item{num.tran}{numeric, the number to transform the OS.time}

\item{break.x}{numeric, the break x}

\item{main.text}{character, the main text of the plot}

\item{color}{character, the color of the plot}
}
\value{
a list
}
\description{
Survival analysis and plot using existed group information
}
\examples{
\dontrun{
df
     OS.time OS  EGFR  Group 
P1     3450  0    4.6   High   
P2     1280  1    1.0    Low   
P3      455  1    3.0    Low   
res <- Surv_diffG(df, col_name, num.tran=365, break.x=2.5, main.text,
                 color=c("#eaa8a8", "#729fc9"))
res <- Surv_diffG(df, col_name, num.tran=365, break.x=2.5, main.text,
                 color=c("#eaa8a8", "#729fc9", "#f6f1af"))
}
}
