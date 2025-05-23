% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_missing_loci.R
\name{get_missing_loci}
\alias{get_missing_loci}
\title{Identify missing loci by analyzing upstream and downstream gap lengths
if they exceed a specified threshold.}
\usage{
get_missing_loci(rmsk_range, sister_clades, folder, chrom = NULL)
}
\arguments{
\item{rmsk_range}{A GRanges object representing repeatMasker regions.}

\item{sister_clades}{A list of sister-clades.}

\item{folder}{A character string specifying the path to store output files.}

\item{chrom}{A character value specifying the chromosome of interest (default:NULL).}
}
\value{
Saves the computed gap lengths in `large_gap_len` and `buffered_large_gap_len directories` in the `folder` directory.
}
\description{
This function determines missing loci by
computing the gap length for these missing loci using immediate and
buffered (±100bp) flanking regions to reduce potential alignment artifacts.
}
\details{
The buffered gap length calculation prevents biases introduced by
artifacts at the very edge of the reference sequence by considering a
region 100bp upstream/downstream.
}
