% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsea.R
\name{cb_fgsea}
\alias{cb_fgsea}
\title{Perform a scored gene set enrichment analysis}
\usage{
cb_fgsea(
  genes,
  scores,
  genesets,
  min_size = 9,
  max_size = 500,
  n_leading_edge = 10,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{genes}{Character vector of genes that have been scored}

\item{scores}{The scores associated with the genes}

\item{genesets}{List of character vectors of pre-defined gene sets (e.g. pathways)}

\item{min_size}{Minimum size of a gene set to be tested (after intersection with input genes) (default is 9)}

\item{max_size}{Maximum size of a gene set to be tested (after intersection with input genes) (default is 500)}

\item{n_leading_edge}{Number of genes to show in the leadingEdge output column - these are the genes driving the enrichment (default is 10)}

\item{verbose}{Boolean indicating whether to show messages (default is TRUE)}

\item{...}{Optional arguments passed on to \code{\link[fgsea:fgsea]{fgsea::fgsea}}}
}
\value{
A table with GSEA results. Each row corresponds to a tested pathway
}
\description{
Perform a scored gene set enrichment analysis
}
\section{Details}{

This is a light wrapper around \code{\link[fgsea:fgsea]{fgsea::fgsea}}. This functions filters the gene sets
given the background and makes sure that there is empty output if there is nothing to be tested.
Output is ordered based on p-value.

This function will always output a data frame. In the case of no results it will have zero rows.
This behavior makes it compatible with a group_by, summarize dplyr workflow.

To load local copies of pathways from EnrichR, check out the helper
function \code{\link{cb_enrichr_gsets}}.
To load them from the web check out \code{\link[hypeR:enrichr_gsets]{hypeR::enrichr_gsets}}.
}

\section{TODO}{

Add documentation and more options to sanitize the geneset labels.
}

\examples{
\donttest{
genesets <- split(LETTERS, f = round((1:length(LETTERS)) / 4))
names(genesets) <- paste0('gene_set_', 1:length(genesets))
genes <- LETTERS
scores <- 1:(length(LETTERS))

fgsea_res <- cb_fgsea(genes = genes,
                      scores = scores,
                      genesets = genesets,
                      min_size = 2)
}

}
