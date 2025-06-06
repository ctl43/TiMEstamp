% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{keep_b4_the_idx}
\alias{keep_b4_the_idx}
\title{Generate Logical Mask for Values Before an Index}
\usage{
keep_b4_the_idx(idx_list, n)
}
\arguments{
\item{idx_list}{An IntegerList specifying the index positions where FALSE values begin.}

\item{n}{An integer or numeric vector indicating the total number of elements in each list component.}
}
\value{
A LogicalList where TRUE represents values to keep and FALSE represents values to discard.
}
\description{
Creates a logical list indicating which values should be retained (TRUE) before a given index.
}
