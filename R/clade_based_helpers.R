#' Aggregate missing fractions by clade
#'
#' Internal helper that computes clade-wise row means from per-species
#' missing fractions.
#'
#' @param x Numeric matrix of missing fractions (loci × species).
#' @param clades List of character vectors; species per clade.
#' @return Numeric matrix (loci × clades) of clade-wise means.
#' @keywords internal
get_missing_by_clades <- function(x, clades){
  x_by_clade <- sapply(clades, function(y){
    rowMeans(x[,y, drop = FALSE], na.rm = TRUE)
  })
  x_by_clade
}

#' Classify existence by clade
#'
#' Internal helper that aggregates missing fractions by clade and
#' assigns clade-level existence calls with propagation once absence appears.
#'
#' @param x Numeric matrix of missing fractions (loci × species).
#' @param clades List of character vectors; species per clade.
#' @param is_abs_target_loci Logical matrix like \code{x}; mask set to \code{NA} before aggregation.
#' @param cutoff Numeric length-2 increasing vector; thresholds for present/absent/uncertain.
#' @param min_species_considered Integer; minimum species per clade to allow absence propagation.
#' @return Integer matrix (loci × clades) with values 1 present, 0 absent, -1 uncertain.
#' @keywords internal
classify_existence_by_clade <- function(x, clades, is_abs_target_loci, cutoff = c(0.25, 0.75), min_species_considered){
  x[is_abs_target_loci] <- NA
  x_by_clade <- get_missing_by_clades(x, clades)
  x_existence_by_clade <- x_by_clade
  x_existence_by_clade[x_by_clade > cutoff[2]] <- 0
  x_existence_by_clade[x_by_clade <= cutoff[1]] <- 1L
  x_existence_by_clade[cutoff[2] >= x_by_clade & x_by_clade > cutoff[1]] <- -1L
  idx_enough_sp <- which(lengths(clades) >= min_species_considered)[1]
  x_existence_by_clade <- t(apply(x_existence_by_clade, 1, function(y){
    which_missing <- which(y == 0)
    which_missing <- which_missing[which_missing >= idx_enough_sp] # because after clade 6 has more samples 
    if(length(which_missing) > 0){
      y[which_missing[1]:length(y)] <- 0L
    }
    return(y)
  }))
  return(x_existence_by_clade)
}

#' Expand clade-level calls to species order
#'
#' Internal helper that repeats clade columns to species membership and
#' reorders columns to a target species order.
#'
#' @param cladewise Integer matrix (loci × clades) of existence calls.
#' @param sister_clades List of character vectors; species per clade.
#' @param colnames_order Character vector; desired final species column order.
#' @return Integer matrix expanded to species and ordered as requested.
#' @keywords internal
expand_cladewise_to_species <- function(cladewise, sister_clades, colnames_order){
  expanded <- mapply(function(x,y){
    do.call(cbind, replicate(y, x, simplify = FALSE))
  },
  x = asplit(cladewise, 2), y = lengths(sister_clades))
  expanded <- do.call(cbind, expanded)
  colnames(expanded) <- unlist(sister_clades)
  expanded <- expanded[, colnames_order]
  return(expanded)
}

