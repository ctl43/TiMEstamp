#' @title Extract First Element from a List
#' @description Extracts the first element of each sub-list within a list.
#' @param x A list of elements.
#' @param invert If TRUE, returns all but the first element.
#' @return A list with the first element extracted (or removed if invert = TRUE).
first_element <- function(x, invert = FALSE){
  if(length(x) == 0){
    return(x)
  }
  grp <- factor(rep(names(x), lengths(x)), levels = names(x))
  flat <- unlist(x, use.names = FALSE)
  is_dup <- duplicated(grp)
  if(invert){
    return(split(flat[is_dup], grp[is_dup]))
  }else{
    return(split(flat[!is_dup], grp[!is_dup]))
  }
}

#' @title Extract Last Element from a List
#' @description Extracts the last element of each sub-list within a list.
#' @param x A list of elements.
#' @param invert If TRUE, returns all but the last element.
#' @return A list with the last element extracted (or removed if invert = TRUE).
last_element <- function(x, invert = FALSE){
  if(length(x) == 0){
    return(x)
  }
  grp <- factor(rep(names(x), lengths(x)), levels = names(x))
  flat <- unlist(x, use.names = FALSE)
  rev_grp <- rev(grp)
  rev_flat <- rev(flat)
  is_dup <- duplicated(rev_grp)
  if(invert){
    return(split(rev(rev_flat[is_dup]), rev(rev_grp[is_dup])))
  }else{
    return(split(rev(rev_flat[!is_dup]), rev(rev_grp[!is_dup])))
  }
}

#' @title Clean Sequence Levels in a GRanges Object
#' @description Removes unused sequence levels from a GRanges object.
#' @param x A GRanges object.
#' @return A GRanges object with unused sequence levels removed.
#' @importFrom GenomeInfoDb seqlevelsInUse
clean_seqlevels <- function(x){
  seqlevels(x) <- seqlevelsInUse(x)
  return(x)
}

#' @title Convert a Matrix to a List
#' @description Converts a matrix into a list, by rows or by columns.
#' @param x A matrix.
#' @param by_col If TRUE, converts by columns; otherwise, converts by rows.
#' @return A list representation of the matrix.
#' @importFrom IRanges NumericList CharacterList IntegerList LogicalList
mat2list <- function(x, by_col = FALSE){
  what <- class(x[1])
  if(by_col){
    out <- apply(x, 2, function(x)x, simplify = FALSE)
  }else{
    out <- apply(x, 1, function(x)x, simplify = FALSE)
  }

  if(what == "numeric"){
    list_fun <- "NumericList"
  }
  if(what == "character"){
    list_fun <- "CharacterList"
  }
  if(what == "integer"){
    list_fun <- "IntegerList"
  }
  if(what == "logical"){
    list_fun <- "LogicalList"
  }
  get(list_fun)(out)
}

#' @title Remove Metadata from a GRanges Object
#' @description Removes metadata columns from a GRanges object.
#' @param x A GRanges object.
#' @return A GRanges object with metadata removed.
empty_metadata <- function(x){
  mcols(x) <- NULL
  return(x)
}

#' @title Generate Logical Mask for Values Before an Index
#' @description Creates a logical list indicating which values should be retained (TRUE) before a given index.
#' @param idx_list An IntegerList specifying the index positions where FALSE values begin.
#' @param n An integer or numeric vector indicating the total number of elements in each list component.
#' @return A LogicalList where TRUE represents values to keep and FALSE represents values to discard.
#' @importFrom IRanges LogicalList
keep_b4_the_idx <- function(idx_list, n){
  # Purpose:
  # Given a list of indices, this function generates a logical mask where:
  # - TRUE values indicate elements that should be retained.
  # - FALSE values indicate elements that should be ignored after the specified index.

  # Step 1: Determine lengths of TRUE regions (before index)
  true2keep_len <- idx_list - 1
  true2keep_len <- true2keep_len[lengths(true2keep_len) > 0]  # Remove empty cases
  true2keep_grp <- factor(rep(names(true2keep_len), unlist(true2keep_len)), levels = names(idx_list))
  true2keep <- rep(TRUE, sum(unlist(true2keep_len)))  # Generate TRUE mask

  # Step 2: Determine lengths of FALSE regions (after index)
  false2append_len <- n - idx_list + 1L
  false2append_len <- false2append_len[lengths(false2append_len) > 0]  # Remove empty cases
  false2append_grp <- factor(rep(names(false2append_len), unlist(false2append_len)), levels = names(idx_list))
  false2append <- rep(FALSE, sum(unlist(false2append_len)))  # Generate FALSE mask

  # Step 3: Combine TRUE and FALSE values into a LogicalList
  is_info2keep <- LogicalList(split(c(true2keep, false2append), c(true2keep_grp, false2append_grp)))

  return(is_info2keep)
}

#' @title Assign FALSE Values After First Encountered FALSE
#' @description Modifies a logical list such that all values after the first encountered FALSE are also set to FALSE.
#' @param pass_list A LogicalList where TRUE represents valid entries and FALSE represents invalid entries.
#' @param to_be_ignored A LogicalList indicating positions that should be ignored when determining the first FALSE.
#' @return A LogicalList where elements after the first detected FALSE remain FALSE.
#' @importFrom IRanges IntegerList LogicalList
assign_falses_after_the_first <- function(pass_list, to_be_ignored){
  # Purpose:
  # This function ensures that once a FALSE is encountered in a vector, all subsequent elements are also set to FALSE.
  # However, if a position is marked in `to_be_ignored`, the first FALSE is ignored, and the function behaves as if
  # the first occurrence was still TRUE.

  # Example 1: Standard Case (no ignored values)
  # Input:   T  T  T  F  T  F  F  T
  # Output:  T  T  T  F  F  F  F  F

  # Example 2: Case where the first value is FALSE
  # Input:   F  T  T  F  T  F  F  T
  # Output:  F  F  F  F  F  F  F  F

  # Example 3: Handling `to_be_ignored`
  # `to_be_ignored`:  T  F  F  F  F  F  F  F (First case is ignored even if FALSE)
  # Input:            F  T  T  F  T  F  F  T
  # Output:           F  T  T  F  F  F  F  F

  # Identify the first occurrence of FALSE that is not ignored
  first_false_idx <- IntegerList(first_element(which(!pass_list & !to_be_ignored)))

  # Generate a logical list to determine values to keep
  is_info2keep <- keep_b4_the_idx(first_false_idx, lengths(pass_list))

  # Identify the length of TRUE values to retain before the first FALSE
  true2keep_len <- first_false_idx - 1
  true2keep_len <- true2keep_len[lengths(true2keep_len) > 0]
  true2keep_grp <- factor(rep(names(true2keep_len), unlist(true2keep_len)), levels = names(first_false_idx))
  true2keep <- rep(TRUE, sum(unlist(true2keep_len)))

  # Identify the length of FALSE values to be assigned after the first FALSE
  false2append_len <- lengths(pass_list) - first_false_idx + 1L
  false2append_len <- false2append_len[lengths(false2append_len) > 0]
  false2append_grp <- factor(rep(names(false2append_len), unlist(false2append_len)), levels = names(first_false_idx))
  false2append <- rep(FALSE, sum(unlist(false2append_len)))

  # Merge TRUE and FALSE values into a logical list
  tf <- c(true2keep, false2append)
  is_info2keep <- LogicalList(split(tf, c(true2keep_grp, false2append_grp)))

  # Handle cases where no FALSE is found in (!pass_list & !to_be_ignored)
  has_no_false <- all(!(!pass_list & !to_be_ignored))
  is_info2keep[has_no_false] <- pass_list[has_no_false] > -1

  # Extract valid values
  info2keep <- pass_list[is_info2keep]
  info2keep_grp <- factor(rep(names(info2keep), lengths(info2keep)), levels = names(info2keep))
  corrected_info <- unname(c(unlist(info2keep, use.names = FALSE), false2append))

  # Reconstruct the corrected logical list
  corrected_false <- split(corrected_info, c(info2keep_grp, false2append_grp))
  corrected_false <- LogicalList(corrected_false)

  # Assign names back to the logical list
  grp <- factor(rep(names(corrected_false), lengths(corrected_false)), levels = names(corrected_false))
  corrected_false <- unlist(corrected_false)
  names(corrected_false) <- names(unlist(pass_list, use.names = FALSE))
  corrected_false <- LogicalList(split(corrected_false, grp))

  return(corrected_false)
}


#' @export
#' @importFrom S4Vectors split
#' @importFrom BiocGenerics unstrand
merge_granges <- function(x, tol, ignore_strand = TRUE){
  .get_list_grp <- function(x, as_factor = TRUE){
    tot <- seq_along(x)
    if(as_factor){
      tot <- factor(tot)
    }
    return(rep(tot, lengths(x)))
  }
  
  
  is_grlist <- class(x)=="CompressedGRangesList"
  if(is_grlist){
    grp <- .get_list_grp(x, as_factor = TRUE)
    gr <- unlist(x)
  }else{
    gr <- x
    grp <- rep(1, length(gr), as_factor = TRUE)
  }
  chrom <- as.integer(seqnames(gr))
  start <- as.integer(start(gr))
  end <- as.integer(end(gr))
  
  if(ignore_strand){
    strand <- rep(1, length(gr), as_factor = TRUE)
    gr <- unstrand(gr)
  }else{
    strand <- as.integer(factor(strand(gr), levels = c("+", "-", "*")))
  }
  
  o <- order(grp, strand, chrom, start, end) # MUST BE SORTED LIKE THIS
  idx <- rep(0, length(gr))
  idx[o] <- cxx_merge_ranges(chrom[o], start[o], end[o], grp = grp[o], strand = strand[o], tol = tol)
  y <- IntegerList(split(idx, grp))
  members <- split(gr, unlist(idx))
  megred <- unlist(range(members))
  out_grp <- y - unname(cumsum(c(0, head(lengths(unique(y)), -1)))) # converting the idx
  out_grp <- IntegerList(out_grp)
  
  if(is_grlist){
    merged_grp <- unique(y)
    out_merged <- split(megred, .get_list_grp(merged_grp, as_factor = TRUE))
    names(out_grp) <- names(out_merged) <- names(x)
    return(list(idx = out_grp, regions = out_merged))
  }else{
    return(list(idx = unlist(out_grp, use.names = FALSE), regions = megred))
  }
}

#'@export
#' @importFrom BiocGenerics paste
convert_to_character <- function(df){
  df <- DataFrame(df)
  converted <- lapply(df, function(x){
    # obj_class <- class(x)
    if(methods::is(x, "GRangesList")){
      grp <- PartitioningByEnd(x)
      x <- relist(as.character(unlist(x, use.names = FALSE)), grp)
      return(paste(x, collapse = ","))
    }
    
    if(grepl(class(x), "List")){
      return(paste(x, collapse = ","))
    }
    
    if(methods::is(x, "list")){
      return(paste(CharacterList(x), collapse = ","))
    }
    return(as.character(x))
  })
  DataFrame(do.call(cbind, converted))
}

