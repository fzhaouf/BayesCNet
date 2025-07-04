% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunBayesCNet.R
\name{GetCellTypeSpecificNetwork}
\alias{GetCellTypeSpecificNetwork}
\title{Get cell-type-specific network}
\usage{
GetCellTypeSpecificNetwork(
  object,
  celltype,
  min_importance = NULL,
  top_n = NULL
)
}
\arguments{
\item{object}{A BayesCNet object with completed inference}

\item{celltype}{Cell type name}

\item{min_importance}{Minimum importance score threshold (default: NULL)}

\item{top_n}{Return top N connections (default: NULL, return all)}
}
\value{
Data frame with network connections for specified cell type
}
\description{
Get cell-type-specific network
}
