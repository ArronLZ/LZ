% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setting.R
\name{setproxy}
\alias{setproxy}
\title{set http proxy for R session}
\usage{
setproxy(ngate = "http://127.0.0.1:7890")
}
\arguments{
\item{ngate}{the proxy address, for example, if you use clash, the address is default "http://127.0.0.1:7890"}
}
\value{
no return，use unsetproxy to cancel the proxy
}
\description{
set http proxy for R session
}
\examples{
# setproxy()
# setproxy("http://127.0.0.1:7890")
}
\author{
Jiang
}
