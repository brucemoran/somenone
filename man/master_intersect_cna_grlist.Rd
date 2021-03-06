% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/somenone_facets.R
\name{master_intersect_cna_grlist}
\alias{master_intersect_cna_grlist}
\title{Find intersect from list of CNA GRanges}
\usage{
master_intersect_cna_grlist(gr_list, ps_vec, which_genome)
}
\arguments{
\item{ps_vec}{mcols columns to keep p(er) s(ample; appended with list element's name)}

\item{which_genome}{the genome assembly used, "hg19" or "hg38"}

\item{cna_list}{list of named GRanges objects}
}
\value{
GRanges object of all CNA, intersecting or not
        includes the ranges and each mcol from ps_vec, annotated `ps_vec[x].sample`
}
\description{
Find intersect from list of CNA GRanges
}
