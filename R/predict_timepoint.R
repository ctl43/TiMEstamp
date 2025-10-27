#' @title Predict Timepoint of Missing Data
#' @description Predicts the timepoint at which missing data starts based on a binary matrix.
#' @param selected_mat A numeric matrix where values represent data presence.
#' @param cutoff Numeric threshold (default = 0.5) for determining missing values.
#' @return A numeric vector indicating the predicted timepoints for each row.
#' @importFrom IRanges IntegerList NumericList
#' @importFrom BiocGenerics which
#' @details
#' Groups 1 to 5 contain only a limited number of species, making them more prone to random missing data.
#' If more than two values are missing within groups 2 to 5, all subsequent groups are considered missing
#' after the first observed missing value.
#' @export
predict_timepoint <- function(selected_mat, cutoff = 0.5){
  rownames(selected_mat) <- seq_len(nrow(selected_mat))
  original <- selected_mat

  # Convert selected_mat to a binary missing matrix (TRUE = missing, FALSE = present)
  selected_binary <- !(round(selected_mat, 2) >= cutoff)

  # Identify the last missing position for each row
  tp_candidate <- last_element(which(mat2list(selected_binary)))
  tp_candidate <- IntegerList(tp_candidate)
  tp_candidate[lengths(tp_candidate) == 0] <- rep(IntegerList(NA), sum(lengths(tp_candidate) == 0))
  tp_candidate <- unname(unlist(tp_candidate))

  # Adjust missing data predictions based on group-specific tolerances
  to_be_taken_care <- lapply(4:12, function(x){
    to_check_col <- 2:(x-1)

    if(x == 4){
      tol <- 1
    }

    if(x == 5){
      tol <- 1
    }

    if(x >= 5){
      tol <- 2
    }

    return(tp_candidate == x & rowSums(!selected_binary[,to_check_col, drop = FALSE], na.rm = TRUE) > tol)
  })

  to_be_taken_care <- Reduce("|", to_be_taken_care)
  to_be_taken_care[is.na(to_be_taken_care)] <- FALSE
  selected_binary[is.na(selected_binary)] <- FALSE

  # Assign FALSE after the first missing value occurrence
  y <- mat2list(selected_binary)
  z <- is.na(y) > -1
  z[to_be_taken_care] <- (is.na(y) < -1)[to_be_taken_care]
  final_binary <- assign_falses_after_the_first(y, z)

  # Determine final timepoints
  final_tp <- last_element(which(final_binary))
  final_tp <- IntegerList(final_tp)
  final_tp[lengths(final_tp) == 0] <- rep(IntegerList(c(0)), sum(lengths(final_tp) == 0))
  final_tp <- unname(unlist(final_tp))

  return(unname(final_tp))
}

#' @title Remove Values Before NA Timepoints
#' @description Removes predicted timepoints if they occur before missing data.
#' @param tp A numeric vector of predicted timepoints.
#' @param mat A matrix indicating missing values.
#' @return A numeric vector where invalid timepoints are set to NA.
#' @export
remove_if_before_na <- function(tp, mat){
  list_idx <- cumsum(duplicated(is.na(mat2list(mat)) < -1)) + 1
  val_after_tp <- mat2list(mat)[list_idx == tp + 1L]

  # Handling edge cases when tp = 11 or tp = 12
  val_after_tp[tp == 11] <- NumericList(999)
  val_after_tp[lengths(val_after_tp) == 0] <- NumericList(999)

  val_after_tp <- unlist(val_after_tp)
  val_after_is_na <- is.na(val_after_tp)

  # Remove timepoints if the next value is NA
  tp[val_after_is_na] <- NA
  return(tp)
}

#' @title Check Timepoint Confidence
#' @description Determines if the predicted timepoint is reliable.
#' @param tp A numeric vector of predicted timepoints.
#' @param mat A matrix indicating missing data.
#' @param na_is_confi Logical flag (default: TRUE) indicating whether NA values should be considered confident.
#' @return A logical vector indicating whether each timepoint is confident.
#' @export
timepoint_is_confi <- function(tp, mat, na_is_confi = FALSE){
  list_idx <- cumsum(duplicated(is.na(mat2list(mat)) < -1)) + 1
  val_after_tp <- mat2list(mat)[list_idx == tp + 1L]

  # Handling edge cases for final timepoints
  val_after_tp[tp == 11] <- NumericList(999)
  val_after_tp[lengths(val_after_tp) == 0] <- NumericList(999)

  val_after_tp <- unlist(val_after_tp)

  # Determine if the value after tp is >= 0.75 (clean data)
  val_after_is_clean <- val_after_tp >= 0.75
  val_after_is_clean[is.na(val_after_is_clean)] <- na_is_confi

  return(val_after_is_clean)
}
