#' Extract representative flanking segments for target loci
#'
#' Iterates over selected chromosomes and calls
#' \code{extract_flanking_segment_by_chrom()} to derive representative
#' upstream or downstream flanking segments. The procedure:
#' \enumerate{
#'   \item loads per-locus flanking gap lengths from
#'         \file{processed/flanking_gap_len/CHROM.rds};
#'   \item keeps species that are consistently absent at the locus
#'         (per clade- and species-level calls);
#'   \item summarizes flanking lengths by cluster to obtain a representative flanking gap after masking target regions (tolerating small variation);
#'   \item tiles each flank into segments and computes missing fractions
#'         per species and by clade;
#'   \item flags segments with inconsistent clade vs. species evidence.
#' }
#'
#' @param folder Character scalar. Project folder containing required inputs and
#'   where outputs will be written.
#' @param which_side Character scalar; \code{"upstream"} or \code{"downstream"}
#'   to select the 5' or 3' flank, respectively.
#' @param chrom Character vector or \code{NULL}. Chromosomes to process. If
#'   \code{NULL}, discovered from \file{folder/gap} by stripping the \code{.rds}
#'   suffix.
#'
#' @details
#' \strong{Required inputs in \code{folder}:}
#' \itemize{
#'   \item \file{sister_clades.rds}: list of species per clade.
#'   \item \file{processed/flanking_gap_len/CHROM.rds}: data frame with columns
#'         \code{up}, \code{dn} (from earlier steps).
#'   \item \file{processed/cladewise_existence/CHROM.rds} and
#'         \file{processed/loci_existence/CHROM.rds}: existence calls.
#'   \item \file{processed/reference_anno/CHROM.rds}: chromosome-specific
#'         reference annotation used to anchor flanks.
#'   \item \file{gap/CHROM.rds}: per-species gap \code{GRanges} with a
#'         \code{labels} factor for species.
#' }
#'
#' \strong{Outputs written:}
#' \itemize{
#'   \item \file{processed/flanking_rep_len/CHROM.rds}: representative flank length per locus (with clustering metadata).
#'   \item \file{processed/flanking_seg_grp/CHROM.rds}: factor mapping segment rows to parent loci.
#'   \item \file{processed/cleaned_flanking_seg_missing_by_clade/CHROM.rds}: clade-wise missing fractions for flanking segments after masking absent loci.
#'   \item \file{processed/flanking_seg_is_cladewise_unreliable/CHROM.rds}: logical list flagging segments with clade–species inconsistency.
#' }
#'
#' @return Invisibly returns \code{NULL}. Called for its side effects of creating
#'   the files listed above.
#'
#' @seealso \code{\link{extract_flanking_segment_by_chrom}},
#'   \code{\link{classify_existence_by_clade}},
#'   \code{\link{expand_cladewise_to_species}},
#'   \code{\link{get_missing_by_clades}},
#'   \code{\link{compute_missing_portion}}
#'
#' @export
extract_flanking_segment <- function(folder, which_side = "upstream", chrom){
  clade_file <- file.path(folder, "sister_clades.rds")
  if(!file.exists(clade_file)){
    stop("No clade information found. Please run get_sister() with your phylogenetic tree in Newick (NH) format first.")
  }
  
  if(is.null(chrom)){
    chrom <- sub(".rds$","",dir(file.path(folder, "gap")))
  }
  for(i in chrom){
    extract_flanking_segment_by_chrom(CHROM = i, folder = folder, which_side = which_side)
  }
}


#' Extract flanking segments for one chromosome
#'
#' Internal worker for \code{extract_flanking_segment()}. Chooses upstream or
#' downstream flank lengths, filters to consistently absent species, clusters
#' lengths to obtain a representative flank per locus, tiles flanks into segments,
#' computes segment-level missing fractions, aggregates by clade, and flags
#' clade–species inconsistencies.
#'
#' @param folder Character. Project folder.
#' @param CHROM Character. Chromosome identifier.
#' @param which_side Character; \code{"upstream"} or \code{"downstream"}.
#'
#' @keywords internal
#' @noRd
#' @importFrom IRanges PartitioningByEnd relist start width subsetByOverlaps IRanges
#' @importFrom S4Vectors splitAsList runLength runValue
#' @importFrom GenomicRanges flank resize GRanges
extract_flanking_segment_by_chrom <- function(folder, CHROM, which_side){
  flanking_gap_len <- readRDS(file.path(folder, "processed", "flanking_gap_len", paste0(CHROM, ".rds")))
  sister_clades <- readRDS(file.path(folder, "sister_clades.rds"))
  if(which_side == "upstream"){
    flanking_len <- flanking_gap_len$up
  }
  
  if(which_side == "downstream"){
    flanking_len <- flanking_gap_len$dn
  }
  
  cladewise_existence_file <- file.path(folder, "processed", "cladewise_existence", paste0(CHROM, ".rds"))
  loci_existence <- readRDS(file.path(folder, "processed", "loci_existence", paste0(CHROM, ".rds")))
  cladewise_existence <- readRDS(cladewise_existence_file)
  seq_is_absent <- make_sure_consistently_absent(loci_existence, cladewise_existence, sister_clades = sister_clades)
  
  flanking_len[!seq_is_absent] <-  (-9999)
  flanking_rep_len <- get_rep_length(flanking_len, tol = 20, filter_by = "max_count", min_count = 5L) # 20 because of TSD
  
  if(!dir.exists(file.path(folder, "processed", "flanking_rep_len"))){
    dir.create(file.path(folder, "processed", "flanking_rep_len"))
  }
  saveRDS(flanking_rep_len, file.path(folder, "processed", "flanking_rep_len", paste0(CHROM, ".rds")))
  
  # Segment the flank regions
  ### This section is equivalent to IntegerList(Map(c,first_len,up_diff)) or  IntegerList(mapply(c, SIMPLIFY = FALSE, first_len, up_diff))
  flat <- unlist(flanking_rep_len$rep_len, use.names = FALSE)
  flanking_seg_len <- integer(length(flat))
  pe <- PartitioningByEnd(flanking_rep_len$rep_len)
  starts_all <- start(pe)
  first_pos  <- starts_all[width(pe) > 0L] # Start positions in the flat vector for nonempty groups
  flanking_seg_len[first_pos] <- flat[first_pos]
  flanking_seg_len[-first_pos] <- unlist(diff(flanking_rep_len$rep_len), use.names = FALSE)
  flanking_seg_len <- relist(flanking_seg_len, pe)
  
  # Making segment GRanges
  reference_anno <- readRDS(file.path(folder, "processed", "reference_anno", paste0(CHROM, ".rds")))
  if(which_side == "upstream"){
    selected_info <- "fivep_frag"
    expanded <- rep(resize(mcols(reference_anno)[[selected_info]], width = 1L, fix = "start"), lengths(flanking_rep_len$rep_len))
    flanking_range <- flank(expanded, width = unlist(flanking_rep_len$rep_len))
    flanking_seg <- resize(flanking_range, width = 1L, fix = "start")
    flanking_seg <- flank(flanking_seg, width = unlist(flanking_seg_len), start = FALSE)
  }
  
  if(which_side == "downstream"){
    selected_info <- "threep_frag"
    expanded <- rep(resize(mcols(reference_anno)[[selected_info]], width = 1L, fix = "end"), lengths(flanking_rep_len$rep_len))
    flanking_range <- flank(expanded, width = unlist(flanking_rep_len$rep_len))
    flanking_seg <- resize(flanking_range, width = 1L, fix = "end")
    flanking_seg <- flank(flanking_seg, width = unlist(flanking_seg_len), start = TRUE)
  }
  
  # Calculating the missing portion of the segments
  missings <- readRDS(file.path(folder, "gap", paste0(CHROM,".rds")))
  included_sp <- unlist(sister_clades, use.names = FALSE)
  missings <- missings[missings$labels %in% included_sp] # Include the selected species only
  missings$labels <- factor(missings$labels, levels = levels(missings$labels)[levels(missings$labels) %in% included_sp])
  missings <- subsetByOverlaps(missings, flanking_seg + 1000L)
  flanking_seg_missing_portion <- compute_missing_portion(missings, flanking_seg)
  
  is_abs_target_loci <- loci_existence == -2
  
  flanking_seg_grp <- factor(rep(as.character(reference_anno), lengths(flanking_rep_len$rep_len)), levels = as.character(reference_anno))
  if(!dir.exists(file.path(folder, "processed", "flanking_seg_grp"))){
    dir.create(file.path(folder, "processed", "flanking_seg_grp"))
  }
  saveRDS(flanking_seg_grp, file.path(folder, "processed", "flanking_seg_grp", paste0(CHROM, ".rds")))
  flanking_seg_existence_cladewise <- classify_existence_by_clade(flanking_seg_missing_portion, clades = sister_clades, 
                                                                  is_abs_target_loci = is_abs_target_loci[flanking_seg_grp, ], 
                                                                  min_species_considered = 15L)
  
  # For chimera prediction
  flanking_seg_missing_portion[is_abs_target_loci[flanking_seg_grp, ]] <- NA
  flanking_seg_missing_by_clade <- get_missing_by_clades(flanking_seg_missing_portion, clades = sister_clades)
  if(!dir.exists(file.path(folder, "processed", "cleaned_flanking_seg_missing_by_clade"))){
    dir.create(file.path(folder, "processed", "cleaned_flanking_seg_missing_by_clade"))
  }
  saveRDS(flanking_seg_missing_by_clade, file.path(folder, "processed", "cleaned_flanking_seg_missing_by_clade", paste0(CHROM, ".rds")))
  
  # Flagging those alignment at the flanking region that is not reliable
  expanded_flanking_cladewise_existence <- expand_cladewise_to_species(flanking_seg_existence_cladewise, sister_clades = sister_clades, colnames_order = colnames(flanking_seg_missing_portion))
  flanking_seg_not_match_count <- rowSums(flanking_seg_missing_portion < 0.75 & expanded_flanking_cladewise_existence == 0, na.rm = TRUE)
  flanking_seg_not_match_cutoff <- ceiling(rowSums(expanded_flanking_cladewise_existence == 0, na.rm = TRUE) * 0.1)
  flanking_seg_is_cladewise_unreliable <- flanking_seg_not_match_count >= flanking_seg_not_match_cutoff
  flanking_seg_is_cladewise_unreliable <- splitAsList(flanking_seg_is_cladewise_unreliable, flanking_seg_grp)
  
  if(!dir.exists(file.path(folder, "processed", "flanking_seg_is_cladewise_unreliable"))){
    dir.create(file.path(folder, "processed", "flanking_seg_is_cladewise_unreliable"))
  }
  saveRDS(flanking_seg_is_cladewise_unreliable, file.path(folder, "processed", "flanking_seg_is_cladewise_unreliable", paste0(CHROM, ".rds")))
}

#' Identify cells absent in both species- and clade-level calls
#'
#' Internal helper that returns a logical matrix where both species-level
#' (\code{x == 0}) and clade-level (\code{evowise == 0}) calls indicate absence.
#'
#' @param x Integer matrix; species-level calls (1 present, 0 absent, -1 uncertain).
#' @param evowise Integer matrix; clade-level calls aligned with \code{sister_clades}.
#' @param sister_clades List of character vectors; species per clade.
#' @return Logical matrix aligned to \code{x}.
#' @keywords internal
#' @noRd
make_sure_consistently_absent <- function(x, evowise, sister_clades){
  evowise[is.na(evowise)] <- 9999 #just to mask NA, so the matrix is filterable.
  evowise <- asplit(evowise, 2)
  is_abs_in_both <- mapply(function(y, z){
    q <- x[, y, drop = FALSE]
    p <- sweep(q==0L, 1, z == 0L, FUN = '&')
    return(p)
  }, y = sister_clades, z = evowise, SIMPLIFY = FALSE)

  is_abs_in_both <- do.call(cbind, is_abs_in_both)
  is_abs_in_both <- is_abs_in_both[, colnames(x)]
  return(is_abs_in_both)
}


#' Representative flank length by clustering
#'
#' Internal helper that clusters per-species flank lengths for each locus with
#' a tolerance \code{tol}, then selects a representative length per locus
#' using \code{filter_by} and \code{min_count}. Returns the selected lengths
#' together with clustering metadata.
#'
#' @param x Numeric matrix (loci × species) of flank lengths.
#' @param tol Integer scalar; clustering tolerance (bp).
#' @param filter_by Character; one of \code{"max_count"} (default) or another
#'   field in the returned object to filter on.
#' @param min_count Integer; minimum cluster count for a locus to be retained.
#'
#' @return A \code{List} with elements:
#'   \itemize{
#'     \item \code{rep_len}: \code{IntegerList} of representative lengths per locus.
#'     \item \code{group_count}: \code{RleList} of cluster counts per locus.
#'     \item \code{max_count}: \code{IntegerList} of maximum counts per locus.
#'     \item \code{selected_cluster}: \code{RleList} of selected cluster indices.
#'     \item \code{cluster_grp}: integer matrix mapping species to clusters.
#'   }
#'
#' @keywords internal
#' @noRd
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges RleList IntegerList
#' @importFrom S4Vectors runLength runValue
get_rep_length <- function(x, tol = 25, filter_by = "max_count", min_count = 5L){
  # Get the representative flanking lengths of each cluster
  dim_names <- dimnames(x)
  x <- mat2list(x)
  
  # Cluster lengths
  grp <- factor(rep(seq_along(x), lengths(x)), levels = seq_along(x))
  flat <- unlist(x, use.names = FALSE)
  temp_gr <- GRanges("x", IRanges(flat, flat))
  merged <- merge_granges(split(temp_gr, grp), tol = tol)
  cluster_grp <- do.call(rbind, merged$idx)
  dimnames(cluster_grp) <- dim_names
  
  # Count the number in each cluster
  idx_count <- RleList(sort(merged$idx))
  names(idx_count) <- names(x)
  
  # Forming sub group
  total_idx <- runValue(idx_count)
  sub_grp_levels <- paste0(rep(names(x), lengths(total_idx)), "_", unlist(total_idx))
  sub_grp <- factor(paste0(rep(names(x), lengths(x)), "_", unlist(merged$idx)), levels = sub_grp_levels)
  
  # Count the occurrence of length in each sub group 
  final_rle <- RleList(sort(IntegerList(split(unname(unlist(x, use.names = FALSE)), sub_grp))))
  final_count <- runLength(final_rle)
  rep_len <- unlist(first_element(runValue(final_rle)[final_count == max(final_count)]))
  max_count <- unlist(first_element(final_count[final_count == max(final_count)]))
  
  # Regroup preventative lengths
  rep_len <- IntegerList(split(unname(rep_len), factor(sub("_.*", "", names(rep_len)), levels = names(x))))
  max_count <- IntegerList(split(unname(max_count), factor(sub("_.*", "", names(max_count)), levels = names(x))))
  group_count <- runLength(idx_count)
  
  # Tidy up the result
  result <- List(rep_len = rep_len, group_count = group_count, max_count = max_count, selected_cluster = total_idx, cluster_grp = cluster_grp)
  pass <- result[[filter_by]] >= min_count & rep_len >= 0
  result$rep_len <- result$rep_len[pass]
  result$group_count <- result$group_count[pass]
  result$max_count <- result$max_count[pass]
  result$selected_cluster <- result$selected_cluster[pass]
  
  return(result)
}
