% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_alignment_coverage.R
\name{get_alignment_coverage}
\alias{get_alignment_coverage}
\title{Computing alignment coverage by clades}
\usage{
get_alignment_coverage(
  tree_file,
  rmsk_file,
  folder,
  consider_missing_loci = TRUE,
  chrom = NULL
)
}
\arguments{
\item{tree_file}{A character string specifying the path to the phylogenetic tree file.}

\item{rmsk_file}{A character string specifying the path to the repeat-masked annotation file.}

\item{folder}{A character string specifying the directory where output files will be stored.}

\item{consider_missing_loci}{A logical value indicating whether missing loci should be considered (default: TRUE).}

\item{chrom}{A character value specifying the chromosome of interest (default:NULL).}
}
\value{
No explicit return value; results are stored in the specified folder.
}
\description{
This function calculates alignment coverage using repeat-masked (rmsk) regions
and phylogenetic tree information. It determines the sister-clades level alignment
coverage of the provided repeat annotation. Additionally, it will compute
if the region is a missing loci if specified.
}
