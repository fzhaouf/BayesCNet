% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{process_gene_inference}
\alias{process_gene_inference}
\title{Process a single gene for Bayesian inference}
\usage{
process_gene_inference(
  gene_idx,
  gene_info,
  data_rna,
  data_atac,
  cell_types_vector,
  Sigma2_est,
  window = 250000,
  regularization = 1e-06
)
}
\arguments{
\item{gene_idx}{Gene index}

\item{gene_info}{Data frame with gene information}

\item{data_rna}{RNA data matrix}

\item{data_atac}{ATAC data matrix}

\item{cell_types_vector}{Cell type assignments}

\item{Sigma2_est}{Estimated covariance matrix}

\item{window}{Regulatory window size}

\item{regularization}{Regularization parameter}
}
\value{
Data frame with inference results
}
\description{
Process a single gene for Bayesian inference
}
\keyword{internal}
