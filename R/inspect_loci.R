#' Inspect missing loci across multiple chromosomes
#'
#' This function is a wrapper that iterates over chromosomes in a project
#' folder and calls [inspect_loci_by_chrom()] for each one. It checks for
#' the presence of sister-clade information, discovers chromosomes to process
#' if none are explicitly provided, and passes parameters downstream.
#'
#' @param folder Character string. Path to the project folder containing
#'   required subdirectories such as \code{gap/} and \code{processed/}.
#' @param selected_info Character or \code{NULL}. If not \code{NULL}, specifies
#'   which metadata column of \code{reference_anno} to use as annotation.
#'   If \code{NULL}, the \code{GRanges/GRangesList} object itself is used.
#' @param min_gap_len_for_abs Integer scalar (default \code{1000L}).
#'   Minimum flanking gap length (in base pairs) used to classify loci as absent.
#' @param chrom Character vector or \code{NULL}. Names of chromosomes to process.
#'   If \code{NULL}, chromosomes are automatically discovered from filenames
#'   under \code{folder/gap}.
#' @param which_side Character string specifying which flanking side(s) 
#'   should determine absence calls for the locus. 
#'   Must be one of \code{"upstream"}, \code{"downstream"}, or \code{"both"}.
#'   - \code{"upstream"}: absence is called if the upstream flank is missing.  
#'   - \code{"downstream"}: absence is called if the downstream flank is missing.  
#'   - \code{"both"}: absence is called only if both flanks are missing.
#' @param buffered_size Integer vector (default \code{c(100, 250)}). Additional
#'   buffer sizes (in base pairs) to use when computing buffered flanking gap lengths.
#'
#' @return Invisibly returns \code{NULL}. As a side effect, calls
#'   [inspect_loci_by_chrom()] for each chromosome, which saves
#'   multiple RDS files to the \code{processed/} subdirectory of \code{folder}.
#'
#' @details
#' Inputs expected in \file{folder}:
#' \itemize{
#'   \item \file{processed/reference_anno/CHROM.rds}: per-locus reference annotation (GRanges/GRangesList).
#'   \item \file{sister_clades.rds}: named list of character vectors defining sister clades.
#'   \item \file{gap/CHROM.rds}: per-species gap GRanges for missing alignment segments.
#'   \item \file{processed/missing_portion/CHROM.rds}: numeric matrix of missing fractions per locus by species.
#' }
#'
#' Files written to \file{folder/processed/}:
#' \itemize{
#'   \item \file{flanking_gap_len/CHROM.rds}: data.frame of unbuffered flanking gap lengths (\code{up}, \code{dn}).
#'   \item \file{cladewise_existence/CHROM.rds}: integer matrix of clade-level existence calls (1 present, 0 absent, -1 uncertain).
#'   \item \file{loci_existence/CHROM.rds}: integer matrix of species-level existence after applying clade masks; \code{-2} marks loci missing the whole segment by flank rule.
#'   \item \file{cleaned_missing_portion/CHROM.rds}: missing-fraction matrix with \code{NA} at loci called absent by flank rule.
#'   \item \file{cleaned_missing_by_clade/CHROM.rds}: row means of cleaned missing fractions aggregated by clade.
#'   \item \file{is_cladewise_unreliable/CHROM.rds}: logical vector flagging loci with inconsistency between clade-level absence and species-level evidence.
#' }
#'
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols isSorted
#' @importFrom XVector compact
#' @export
inspect_loci <- function(folder, selected_info = "fivep_frag", min_gap_len_for_abs = 1000L, chrom = NULL, which_side = "upstream", buffered_size = c(100, 250)){
  clade_file <- file.path(folder, "sister_clades.rds")
  if(!file.exists(clade_file)){
    stop("No clade information found. Please run get_sister() with your phylogenetic tree in Newick (NH) format first.")
  }
  
  if(is.null(chrom)){
    chrom <- sub(".rds$","",dir(file.path(folder, "gap")))
  }
  for(i in chrom){
    inspect_loci_by_chrom(CHROM = i, folder = folder, selected_info = selected_info, min_gap_len_for_abs = min_gap_len_for_abs, buffered_size = buffered_size, which_side = which_side)
  }
}

#' Inspect one chromosome and write processed artifacts
#'
#' Internal worker for \code{inspect_loci()}: computes flanking-gap lengths,
#' reconciles species- and clade-level existence calls, and saves per-step RDS files.
#'
#' @param CHROM Character; chromosome identifier.
#' @param folder Character; project root directory.
#' @param min_gap_len_for_abs Integer; flank-length threshold (bp).
#' @param selected_info Character or NULL; reference annotation selector.
#' @param buffered_size Integer vector; extra symmetric flank buffers (bp).
#' @param which_side Character; one of "upstream", "downstream", "both".
#' @param min_available_data Numeric; reserved for future use.
#'
#' @keywords internal
#' @noRd
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols isSorted
#' @importFrom XVector compact
inspect_loci_by_chrom <- function(CHROM, folder, min_gap_len_for_abs, selected_info, buffered_size, which_side = "up", min_available_data = 0.1){
  if(!dir.exists(folder)){
    stop("Output folder does not exist.")
  }
  reference_anno <- readRDS(file.path(folder, "processed", "reference_anno", paste0(CHROM, ".rds")))
  sister_clades <- readRDS(file.path(folder, "sister_clades.rds"))
  
  # Select useful missing data only
  if(is.null(selected_info)){
    anno_names <- names(reference_anno)
    message("No metadata selected; using the GRanges/GRangesList object provided in reference_file as input.")
  }else{
    anno_names <- names(reference_anno)
    reference_anno <- mcols(reference_anno)[[selected_info]]
    names(reference_anno) <- anno_names
  }
  
  missings <- readRDS(file.path(folder, "gap", paste0(CHROM,".rds")))
  
  if(selected_info == "members"){
    missings <- subsetByOverlaps(missings, range(reference_anno) + min_gap_len_for_abs + max(buffered_size) + 1000L)
  }else{
    missings <- subsetByOverlaps(missings, reference_anno + min_gap_len_for_abs + max(buffered_size) + 1000L)
  }
  
  missings <- XVector::compact(missings)
  
  # Include the selected species only
  included_sp <- unlist(sister_clades, use.names = FALSE)
  missings <- missings[missings$labels %in% included_sp]
  missings$labels <- factor(missings$labels, levels = levels(missings$labels)[levels(missings$labels) %in% included_sp])
  
  # Calculating flanking gap length
  flanking_gap_len <- list()
  flanking_gap_len[[1L]] <- compute_flanking_gap_length(missings, ref = reference_anno)
  
  # Saving this information for getting representative flanking gap length to find chimeric insertions
  if(!dir.exists(file.path(folder, "processed","flanking_gap_len"))){
    dir.create(file.path(folder, "processed", "flanking_gap_len"))
  }
  saveRDS(flanking_gap_len[[1L]], file.path(folder, "processed", "flanking_gap_len", paste0(CHROM, ".rds")))
  
  if(length(buffered_size) != 0){
    for (i in seq_along(buffered_size)){
      flanking_gap_len[[i + 1L]] <- compute_flanking_gap_length(missings, ref = reference_anno + buffered_size[i])
    }
  }

  missing_portion <- readRDS(file.path(folder, "processed", "missing_portion", paste0(CHROM, ".rds")))
  existence <- classify_existence(missing_portion, cutoff = c(0.25, 0.75)) # 1: L1 present 0: L1 absent -1: not sure
  
  # Finding the species that miss the corresponding segment entirely
  if(which_side == "both"){
    is_abs_target_loci_draft <- existence != 1 & Reduce("|", lapply(flanking_gap_len, function(x)(x$up > min_gap_len_for_abs)|(x$dn > min_gap_len_for_abs)))
  }
  
  if(which_side == "upstream"){
    is_abs_target_loci_draft <- existence != 1 & Reduce("|", lapply(flanking_gap_len, function(x)x$up > min_gap_len_for_abs))
  }
  
  if(which_side == "downstream"){
    is_abs_target_loci_draft <- existence != 1 & Reduce("|", lapply(flanking_gap_len, function(x)x$dn > min_gap_len_for_abs))
  }
  
  # Classifying the presence of L1 - by evolution mean
  cladewise_existence <- classify_existence_by_clade(missing_portion, sister_clades, is_abs_target_loci_draft, min_species_considered = 10L)
  if(!dir.exists(file.path(folder, "processed","cladewise_existence"))){
    dir.create(file.path(folder, "processed", "cladewise_existence"))
  }
  saveRDS(cladewise_existence, file.path(folder, "processed", "cladewise_existence", paste0(CHROM, ".rds")))
  
  # Combining the existence matrix computed by missing portion and evolution logic
  existence_combined <- apply_clade_existence_mask(existence, cladewise_existence, sister_clades = sister_clades)
  
  if(which_side == "both"){
    is_abs_target_loci_final <- existence_combined != 1 & Reduce("|", lapply(flanking_gap_len, function(x)(x$up > min_gap_len_for_abs)|(x$dn > min_gap_len_for_abs)))
  }
  
  if(which_side == "upstream"){
    is_abs_target_loci_final <- existence_combined != 1 & Reduce("|", lapply(flanking_gap_len, function(x)x$up > min_gap_len_for_abs))
  }
  
  if(which_side == "downstream"){
    is_abs_target_loci_final <- existence_combined != 1 & Reduce("|", lapply(flanking_gap_len, function(x)x$dn > min_gap_len_for_abs))
  }
  
  # saveRDS(is_abs_target_loci_final, file.path(folder, "processed", CHROM, "is_abs_target_loci_final.rds"))
  
  existence_combined[is_abs_target_loci_final] <- -2 # -2: miss the whole segment
  if(!dir.exists(file.path(folder, "processed","loci_existence"))){
    dir.create(file.path(folder, "processed", "loci_existence"))
  }
  saveRDS(existence_combined, file.path(folder, "processed", "loci_existence", paste0(CHROM, ".rds")))
  
  # Useful data for downstream analysis, for example plotting heatmap or predicting chimera
  cleaned_missing_portion <- missing_portion
  cleaned_missing_portion[is_abs_target_loci_final] <- NA
  if(!dir.exists(file.path(folder, "processed","cleaned_missing_portion"))){
    dir.create(file.path(folder, "processed", "cleaned_missing_portion"))
  }
  saveRDS(cleaned_missing_portion, file.path(folder, "processed", "cleaned_missing_portion", paste0(CHROM, ".rds")))
  
  cleaned_missing_by_clade <- get_missing_by_clades(cleaned_missing_portion, sister_clades)
  if(!dir.exists(file.path(folder, "processed","cleaned_missing_by_clade"))){
    dir.create(file.path(folder, "processed", "cleaned_missing_by_clade"))
  }
  saveRDS(cleaned_missing_by_clade, file.path(folder, "processed", "cleaned_missing_by_clade", paste0(CHROM, ".rds")))
  
  rm(flanking_gap_len); rm(cleaned_missing_portion); rm(missings); rm(cleaned_missing_by_clade); gc()
  
  expanded_cladewise_existence <- expand_cladewise_to_species(cladewise_existence, sister_clades = sister_clades, colnames_order = colnames(existence_combined))
  
  # # By coherence between evolution wise and individual
  expanded_cladewise_existence[is_abs_target_loci_final] <- NA
  is_cladewise_absent <- expanded_cladewise_existence == 0
  not_match_count <- rowSums(missing_portion < 0.75 & is_cladewise_absent, na.rm = TRUE)
  not_match_cutoff <- ceiling(rowSums(is_cladewise_absent, na.rm = TRUE) * 0.1)
  is_cladewise_unreliable <- not_match_count >= not_match_cutoff
  
  if(!dir.exists(file.path(folder, "processed","is_cladewise_unreliable"))){
    dir.create(file.path(folder, "processed", "is_cladewise_unreliable"))
  }
  saveRDS(is_cladewise_unreliable, file.path(folder, "processed", "is_cladewise_unreliable", paste0(CHROM, ".rds")))
}

#' Apply clade-level absence mask
#'
#' Internal helper to reconcile species- and clade-level existence calls.
#' Sets species calls to -1 where clade-level absence contradicts them.
#' @param x Integer matrix of species-level calls (1 present, 0 absent, -1 uncertain).
#' @param evowise Integer matrix of clade-level calls (1 present, 0 absent, -1 uncertain).
#' @param sister_clades List of character vectors giving species in each clade.
#' @keywords internal
apply_clade_existence_mask <- function(x, evowise, sister_clades){
  # If it is contradict, assign -1 indicate uncertain
  evowise[is.na(evowise)] <- 9999 #just to mask NA, so the matrix is filterable. 
  evowise <- asplit(evowise, 2)
  cleaned <- mapply(function(y, z){
    q <- x[, y, drop = FALSE]
    p <- sweep(q!=0L, 1, z == 0, FUN = '&')
    q[p] <- -1L
    return(q)
  }, y = sister_clades, z = evowise)
  
  cleaned <- do.call(cbind, cleaned)
  cleaned <- cleaned[, colnames(x)]
  return(cleaned)
}


# mask_insufficient_clade_data <- function(missing_by_clade, loci_existence, sister_clades, min_available_data){
#   n_available_data_by_clade <- sapply(sister_clades, function(y)rowSums(!loci_existence[,y, drop = FALSE]))
#   n_data_cutoff <- as.integer(lengths(sister_clades) * min_available_data)
#   n_data_cutoff[n_data_cutoff <= 2]  <- 0 # If clades have less than 20 assemblies, do nothing
#   n_data_cutoff <- t(replicate(nrow(n_available_data_by_clade), n_data_cutoff)) # can be improved
#   missing_by_clade[n_data_in_clades < n_data_cutoff] <- NA
#   return(missing_by_clade)
# }

# apply_clade_mask <- function(missings, missing_by_clade, sister_clades){
#   by_clade <- asplit(missing_by_clade, 2)
#   masked_missings <- lapply(sister_clades, function(x){
#     y <- missings[, x, drop = FALSE]
#     y[is.na(by_clade)] <- NA
#     return(y)
#   })
#   return(do.call(cbind, masked_missings))
# }



