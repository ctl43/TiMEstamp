% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_timepoint.R
\name{predict_timepoint}
\alias{predict_timepoint}
\title{Predict Timepoint of Missing Data}
\usage{
predict_timepoint(selected_mat, cutoff = 0.5)
}
\arguments{
\item{selected_mat}{A numeric matrix where values represent data presence.}

\item{cutoff}{Numeric threshold (default = 0.5) for determining missing values.}
}
\value{
A numeric vector indicating the predicted timepoints for each row.
}
\description{
Predicts the timepoint at which missing data starts based on a binary matrix.
}
\details{
Groups 1 to 5 contain only a limited number of species, making them more prone to random missing data.
If more than two values are missing within groups 2 to 5, all subsequent groups are considered missing
after the first observed missing value.
}
