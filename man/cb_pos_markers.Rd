% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_expression.R
\name{cb_pos_markers}
\alias{cb_pos_markers}
\title{Find positive cluster markers}
\usage{
cb_pos_markers(counts, grouping, FDR_cutoff = 0.05)
}
\arguments{
\item{counts}{Count matrix with features as rows and cells as columns}

\item{grouping}{Factor or character vector indicating the group membership of the cells}

\item{FDR_cutoff}{Remove results with FDR below this value; set to NULL or Inf to keep all
results; default is 0.05}
}
\value{
A data frame of results
}
\description{
Find positive cluster markers
}
\section{Details}{

For each group, this creates pseudobulk data for the cells in that group and
for the cells that are not in the group. The two pseudobulk sets are compared
against each other using edgeR. Only results with positive fold-change are kept
and FDR values are calculated after all groups have been processed.
}

\examples{
\donttest{
# Coming soon
}

}
