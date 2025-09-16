#' Classify species- or locus-level existence from missing fractions
#'
#' Converts missing fraction values into categorical existence calls
#' using lower and upper cutoffs. Values below the lower cutoff are
#' called present, values above the upper cutoff are called absent,
#' and values between cutoffs are called uncertain.
#'
#' @param missing_portion Numeric vector or matrix. Missing fraction(s)
#'   per locus or species. Values are typically proportions in \[0,1].
#' @param cutoff Numeric vector of length two, sorted in increasing order.
#'   Defines thresholds for classification:
#'   \itemize{
#'     \item \code{< cutoff[1]} → present (\code{1L})
#'     \item \code{> cutoff[2]} → absent (\code{0L})
#'     \item \code{between cutoff[1] and cutoff[2]} → uncertain (\code{-1L})
#'   }
#'   Default is \code{c(0.05, 0.95)}.
#'
#' @return An integer vector or matrix with the same dimensions and names
#'   as \code{missing_portion}. Values are:
#'   \itemize{
#'     \item \code{1} = present
#'     \item \code{0} = absent
#'     \item \code{-1} = uncertain
#'   }
#'
#' @importFrom S4Vectors isSorted

classify_existence <- function(missing_portion, cutoff = c(0.05, 0.95)){
  # cutoff here means missing portion
  # 0.05 means missing 5% of the alignemnt
  # 0.95 means missing 95% of the alignment
  if(!isSorted(cutoff)){
    stop("Cutoff must be sorted")
  }
  existence <- missing_portion
  existence[missing_portion < cutoff[1]] <- 1L # exist
  existence[missing_portion > cutoff[2]] <- 0L # not exist
  existence[missing_portion >= cutoff[1] & missing_portion <= cutoff[2]] <- -1L # not sure
  if("matrix" %in% class(existence)){
    dimnames(existence) <- dimnames(missing_portion)
  }
  if(any(c("numeric", "integer") %in% class(existence))){
    names(existence) <- names(missing_portion)
  }
  return(existence)
}
