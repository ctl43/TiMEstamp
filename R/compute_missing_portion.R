#' Overlap-based missing-fraction matrix
#'
#' Internal helper that computes per-species missing fractions by intersecting
#' gap intervals with reference intervals and optionally normalizing by reference length.
#'
#' @param x GRanges of gap intervals with a \code{labels} factor for species.
#' @param ref GRanges or CompressedGRangesList; reference intervals.
#' @param name Optional character vector of row names for the output matrix.
#' @param normalized_by_ref_length Logical; divide by reference length if \code{TRUE}.
#'
#' @return Numeric matrix with rows as loci and columns as species.
#'
#' @family missing-fraction
#' @keywords internal
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges pintersect width IntegerList
compute_missing_portion <- function(x, ref, name = NULL, normalized_by_ref_length = TRUE){
  ref_is_CompressedGRangesList <- FALSE
  if(class(ref) == "CompressedGRangesList"){
    ref_is_CompressedGRangesList <- TRUE
    list_grp <- factor(rep(seq_along(ref), lengths(ref)), levels = seq_along(ref))
    ref_length <- sum(width(ref))
    ref_names <- names(ref)
    ref <- unlist(ref, use.names = FALSE)
  }else{
    ref_names <- as.character(ref)
    ref_length <- width(ref)
  }
  ol <- findOverlaps(x, ref)
  intersected <- pintersect(x[ol@from], ref[ol@to])
  splitted_to <- split(factor(ol@to, levels = seq_along(ref)), x[ol@from]$labels)
  splitted_width <- split(width(intersected), x[ol@from]$labels)
  out <- mapply(function(y,z){
    sum(IntegerList(split(y, z)))
  }, y = splitted_width, z = splitted_to)
  if(!grepl("matrix", class(out)[1])){
    out <- matrix(out, 1, , )
  }
  if(ref_is_CompressedGRangesList){
    out <- rowsum(out, list_grp, reorder = FALSE)
  }
  
  if(normalized_by_ref_length){
    normalized <- apply(out, 2, function(x)x/ref_length)
  }else{
    normalized <- out
  }
  
  if(!grepl("matrix", class(normalized)[1])){
    normalized <- matrix(normalized, 1, , )
  }
  if(is.null(name)){
    rownames(normalized) <- ref_names
  }else{
    rownames(normalized) <- name
  }
  colnames(normalized) <- levels(x$labels)
  return(normalized)
}

