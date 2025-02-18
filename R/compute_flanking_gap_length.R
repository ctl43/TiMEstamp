#' Compute Flanking Gap Lengths of Loci That Do Not Have The Insertion
#'
#' This function calculates the upstream and downstream gap lengths of locs that do not have the TE insertion
#' identifies flanking regions based on overlaps and extracts gaps using a logic similar to `setdiff` in `GenomicRanges`.
#'
#' @param missings A `GRanges` object representing missing regions in the genome.
#' @param ref A `GRanges` object representing the reference TE annotations.
#'
#' @return A list containing:
#'   - `up`: A matrix of upstream gap lengths grouped by species.
#'   - `dn`: A matrix of downstream gap lengths grouped by species.
#'
#' @export
#' @importFrom GenomicRanges findOverlaps GRanges seqnames start end resize flank
#' @importFrom S4Vectors mcols<-
#' @importFrom IRanges IntegerList ranges
#' @importFrom stats setNames
#' @importFrom BiocGenerics width
compute_flanking_gap_length <- function(missings, ref) {
  # Identify overlaps between missing regions and reference TE annotations
  # This helps extract the flanking sequences outside the TE annotations

  ol <- findOverlaps(missings, ref + 1L)
  # ---------------------------- This is equavilent to setdiff
  # set diff, masking the refernce sequneces in missings to get the flanking regions of TE
  # Case 1
  # ----------------   Missing regions
  #               ---- TE annotation
  # --------------     Result
  # Case 2
  #   ---------------- Missing regions
  # ----       TE annotation
  #     -------------- Result
  # Case 3
  #       ----       Missing regions
  # ---------------- TE annotation
  #                  Result
  # Case 4
  # ---------------- Missing regions
  #       ----       TE annotation
  # ------    ------ Result
  # Extract coordinate values
  x_labels <- missings$labels[ol@from]
  x_start <- start(missings)[ol@from]
  x_end <- end(missings)[ol@from]
  y_start <- start(ref)[ol@to]
  y_end <- end(ref)[ol@to]
  x_seqnames <- seqnames(missings)[ol@from]

  # Determine cases based on overlap patterns
  y_start_larger <- x_start < y_start  # Case where the missing region starts before TE
  y_end_smaller <- y_end < x_end        # Case where the TE ends before the missing region

  # Define different overlap cases
  case1 <- y_start_larger & !y_end_smaller # Missing starts before TE but ends at TE
  case2 <- !y_start_larger & y_end_smaller # Missing starts inside TE but ends outside
  case3 <- !y_start_larger & !y_end_smaller # Missing is completely inside TE (ignored)
  case4 <- y_start_larger & y_end_smaller  # Missing extends before and after TE
  # -------------------------------------------------------------------------------------------
  # Process Case 1: Extract the left flanking region of TE
  if (length(case1) > 0) {
    case1_out <- GRanges(paste0(x_seqnames[case1], "_", ol@to[case1]),
                         IRanges(x_start[case1], y_start[case1] - 1L))
  } else {
    case1_out <- GRanges()
  }

  # Process Case 2: Extract the right flanking region of TE
  if (length(case2) > 0) {
    case2_out <- GRanges(paste0(x_seqnames[case2], "_", ol@to[case2]),
                         IRanges(y_end[case2] + 1L, x_end[case2]))
  } else {
    case2_out <- GRanges()
  }

  # Process Case 4: Extract both flanking regions
  if (length(case4) > 0) {
    case4_out_p1 <- GRanges(paste0(x_seqnames[case4], "_", ol@to[case4]),
                            IRanges(x_start[case4], y_start[case4] - 1L))
    case4_out_p2 <- GRanges(paste0(x_seqnames[case4], "_", ol@to[case4]),
                            IRanges(y_end[case4] + 1L, x_end[case4]))
  } else {
    case4_out_p1 <- case4_out_p2 <- GRanges()
  }

  case4_out <- c(case4_out_p1, case4_out_p2)

  # Assign labels to extracted regions
  case1_out$labels <- x_labels[case1]
  case2_out$labels <- x_labels[case2]
  case4_out$labels <- c(x_labels[case4], x_labels[case4])

  # Combine all cases into a final flanking region output
  out <- c(case1_out, case2_out, case4_out)
  grp <- c(ol@to[case1], ol@to[case2], ol@to[case4], ol@to[case4])

  ### Compute Upstream Gaps ###
  ref_up <- ref
  mcols(ref_up) <- NULL  # Remove metadata to simplify processing
  ref_up <- flank(resize(ref_up, width = 1), 10L)  # Extract upstream region (10bp)
  ref_up <- GRanges(paste0(seqnames(ref_up), "_", seq_along(ref_up)), ranges(ref_up))

  # Find overlaps between missing regions and upstream reference sequences
  ol_up <- findOverlaps(out, ref_up)
  masked_up <- out[ol_up@from]
  masked_up <- GRanges(gsub("_[0-9]+", "", seqnames(masked_up)), ranges(masked_up), labels = masked_up$labels)

  up_target_grp <- factor(ol_up@to, levels = seq_along(ref_up))
  up_len_by_sp <- split(width(masked_up), masked_up$labels)
  up_grp_by_sp <- split(up_target_grp, masked_up$labels)

  up_length <- mapply(function(y, z) {
    sum(IntegerList(split(y, z)))
  }, y = up_len_by_sp, z = up_grp_by_sp)

  rownames(up_length) <- as.character(ref)

  ### Compute Downstream Gaps ###
  ref_dn <- ref
  mcols(ref_dn) <- NULL
  ref_dn <- flank(resize(ref_dn, width = 1, fix = "end"), 10L, start = FALSE)  # Extract downstream region (10bp)
  ref_dn <- GRanges(paste0(seqnames(ref_dn), "_", seq_along(ref_dn)), ranges(ref_dn))

  ol_dn <- findOverlaps(out, ref_dn)
  masked_dn <- out[ol_dn@from]
  masked_dn <- GRanges(gsub("_[0-9]+", "", seqnames(masked_dn)), ranges(masked_dn), labels = masked_dn$labels)

  dn_target_grp <- factor(ol_dn@to, levels = seq_along(ref_dn))
  dn_len_by_sp <- split(width(masked_dn), masked_dn$labels)
  dn_grp_by_sp <- split(dn_target_grp, masked_dn$labels)

  dn_length <- mapply(function(y, z) {
    sum(IntegerList(split(y, z)))
  }, y = dn_len_by_sp, z = dn_grp_by_sp)

  rownames(dn_length) <- as.character(ref)

  # Clean memory
  gc()

  # Return a list containing upstream and downstream flanking gap lengths
  return(list(up = up_length, dn = dn_length))
}
