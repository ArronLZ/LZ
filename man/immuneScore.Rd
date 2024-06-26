% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/immune.R
\name{immuneScore}
\alias{immuneScore}
\title{immuneScore}
\usage{
immuneScore(
  exprs,
  method,
  tcga_abbr,
  perm = 1000,
  QN = F,
  absolute = T,
  abs_method = "sig.score",
  gsva_sig_list
)
}
\arguments{
\item{exprs}{matrix, gene expression matrix, NOT log-transformed}

\item{method}{character, one of "mcp_counter", "quantiseq", "epic", "xcell", "estimate", "timer", "cibersort", "ssGSEA".}

\item{tcga_abbr}{character, TCGA abbreviation, only used for method = "timer"}

\item{perm}{number, default 1000(recommend 1000), only used for method = "cibersort"}

\item{QN}{logical, default F, RNAseq data: recommend QN = F, otherwise T, only used for method = "cibersort"}

\item{absolute}{absolute logical, if T, return absolute score, only used for method = "cibersort"}

\item{abs_method}{character, default 'sig.score', can be 'no.sumto1', only used for method = "cibersort"}

\item{gsva_sig_list}{list, default NULL, only used for method = "ssGSEA"}
}
\value{
a data frame
}
\description{
Calculate immune score using RNAseq data or mcroarray data.
}
\details{
The input data is a gene × sample gene expression matrix. \cr
In general values should be \cr
 TPM-normalized \cr
 not log-transformed. \cr
For xCell and MCP-counter this is not so important. \cr
 xCell works on the ranks of the gene expression only and 
 MCP-counter sums up the gene expression values.\cr
Rownames are expected to be HGNC gene symbols. \cr
Instead of a matrix, immunedeconv also supports ExpressionSets\cr
=======================================\cr
recommend mcp_counter when comparing between samples
EPIC, quanTIseq, CIBERSORT abs : **both**\cr
CIBERSORT : **between-cell-type comparisons**\cr
MCP-counter, xCell, TIMER, ConsensusTME, ESTIMATE, ABIS: **between-sample comparisons**
}
\examples{
\dontrun{
exprs <- read.table("xx.txt", sep = "\t", row.names = 1, header = T, check.names = F)
exprs \%>\% lzhead
 #          TCGA-05-4249-01A TCGA-05-4250-01A TCGA-05-4382-01A
 # FGB            84.3223886         17.19138     2495.5995163
 # PITX2         127.4491065        125.84088       49.7737300
 # SERPINA4        1.2873647       3021.55644        0.0000000
 }
}
