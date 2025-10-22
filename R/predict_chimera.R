#' Predict chimeric insertions
#'
#' For each selected chromosome, predicts loci with putative chimeric insertions
#' by comparing the inferred timepoint of the locus to the timepoint
#' inferred from its flanking segments, and by requiring consistent support from
#' segment clusters and clade-level reliability filters.
#'
#' @param folder Character scalar. Project folder containing required inputs and
#'   where outputs will be written under \file{processed/}.
#' @param chrom Character vector or \code{NULL}. Chromosomes to process. If
#'   \code{NULL}, discovered from \file{folder/gap} by stripping the \code{.rds} suffix.
#' @param which_side Character; \code{"upstream"} or \code{"downstream"} to choose
#'   the 5' or 3' flank, respectively. Default \code{"upstream"}.
#' @param min_len Integer (bp). Minimum representative flank length to consider
#'   as a candidate segment. Default \code{25}.
#' @param max_len Integer (bp). Maximum representative flank length to consider.
#'   Default \code{1000}. (Current implementation stops if \code{max_len} is below
#'   \code{1000}.)
#' @param min_cov Numeric in \code{[0,1]}. Minimum clade-level coverage agreement
#'   (converted internally to a missingness cutoff as \code{1 - min_cov}) used by
#'   \code{predict_timepoint()}. Default \code{0.65}.
#' @param min_supporting Integer. Minimum total count of supporting segments
#'   across clusters (after filtering) required to call a locus. Default \code{15L}.
#'
#' @details
#' \strong{Inputs expected (created by earlier steps):}
#' \itemize{
#'   \item \file{processed/cleaned_missing_by_clade/CHROM.rds} — clade-wise missing fractions for loci.
#'   \item \file{processed/flanking_rep_len/CHROM.rds} — representative flank lengths and clustering metadata.
#'   \item \file{processed/cleaned_flanking_seg_missing_by_clade/CHROM.rds} — clade-wise missing fractions for flanking segments.
#'   \item \file{processed/flanking_seg_grp/CHROM.rds} — factor mapping segment rows to parent loci.
#'   \item \file{processed/flanking_seg_is_cladewise_unreliable/CHROM.rds} — segment-level reliability flags.
#'   \item \file{processed/is_cladewise_unreliable/CHROM.rds} — locus-level reliability flags.
#'   \item \file{processed/reference_anno/CHROM.rds} — per-chromosome reference loci.
#' }
#'
#' \strong{Outputs:}
#' \itemize{
#'   \item \file{processed/predicted_chimera_<which_side>/CHROM.rds} — a subset of
#'         \code{reference_anno} with additional fields for predicted flanking
#'         coordinates (\code{$flanking_gcoord}) and combined locus range
#'         (\code{$full_range}) for passed loci.
#' }
#'
#' The function validates basic parameter constraints (e.g., \code{min_len > 10},
#' \code{max_len} not below \code{1000}, and \code{max_len > min_len}).
#'
#' @return Invisibly returns \code{NULL}. Called for its side effects (files written).
#'
#' @seealso \code{\link{extract_flanking_segment}}, \code{\link{predict_timepoint}}
#'
#' @export
predict_chimera <- function(folder, chrom = NULL, which_side = "upstream",  min_len = 25, max_len = 1000, min_cov = 0.65, min_supporting = 15L){
  if(is.null(chrom)){
    chrom <- sub(".rds$","",dir(file.path(folder, "gap")))
  }
  for(i in chrom){
    predict_chimera_by_chrom(folder = folder, CHROM = i, which_side = which_side, min_len = min_len, max_len = max_len, min_cov = min_cov, min_supporting = min_supporting)
  }
  out_folder <- paste0("predicted_chimera", "_", which_side)
  all_out <- dir(file.path(folder, "processed", out_folder), full.names = TRUE)
  out <- lapply(all_out, readRDS)
  out <- Reduce(c, out)
  df_out <- convert_to_character(mcols(out))
  out <- cbind(DataFrame(coordinate = as.character(out)), df_out)
  write.table(as.data.frame(out), file.path(folder, paste0(paste0("predicted_chimera", "_", which_side), ".txt")), row.names = FALSE, quote = FALSE, sep = "\t")
}

#' Predict chimeric insertions for one chromosome
#'
#' Internal worker for \code{predict_chimera()}. Derives timepoints for loci and
#' their flanking segments, enforces length and reliability filters, and writes
#' the subset of reference loci that pass chimera criteria.
#'
#' @param folder Character. Project folder.
#' @param CHROM Character. Chromosome identifier.
#' @param which_side Character; \code{"upstream"} or \code{"downstream"}.
#' @param min_len Integer (bp). Minimum representative flank length.
#' @param max_len Integer (bp). Maximum representative flank length.
#' @param min_cov Numeric in \code{[0,1]}. Converted to \code{1 - min_cov} for \code{predict_timepoint()}.
#' @param min_supporting Integer. Minimum supporting segment count per locus.
#'
#' @keywords internal
#' @noRd
#' @importFrom GenomicRanges flank resize
#' @importFrom S4Vectors splitAsList mcols
predict_chimera_by_chrom <- function(folder, CHROM, which_side, min_len, max_len, min_cov, min_supporting){
  
  if(min_len <= 10L){
    stop("min_len must be larger than 10L")
  }
  
  if(max_len < 1000L){
    stop("max_len must be smaller than 1000L")
  }
  
  if(max_len <= min_len){
    stop("max_len must be larger than min_len")
  }
  
  missing_by_clade <- readRDS(file.path(folder, "processed", "cleaned_missing_by_clade", paste0(CHROM, ".rds")))
  tp <- predict_timepoint(missing_by_clade, cutoff = (1-min_cov))
  names(tp) <- rownames(missing_by_clade)
  
  flanking_seg_missing_by_clade <- readRDS(file.path(folder, "processed", "cleaned_flanking_seg_missing_by_clade", paste0(CHROM, ".rds")))
  flanking_seg_grp <- readRDS(file.path(folder, "processed", "flanking_seg_grp", paste0(CHROM, ".rds")))
  flanking_seg_tp <- predict_timepoint(flanking_seg_missing_by_clade, cutoff = (1-min_cov)) # 1/3
  flanking_seg_tp <- splitAsList(flanking_seg_tp, flanking_seg_grp)
  flanking_rep_len <- readRDS(file.path(folder, "processed", "flanking_rep_len", paste0(CHROM, ".rds")))
  
  flanking_length_pass <- flanking_rep_len$rep_len > min_len & flanking_rep_len$rep_len < max_len
  same_tp <- flanking_seg_tp == tp
  flanking_pass_1 <- same_tp & flanking_length_pass # & co_percentage_pass #& up_seg_edge_pass
  flanking_is_zero <- flanking_rep_len$rep_len < 10L
  flanking_pass_2 <- assign_falses_after_the_first(flanking_pass_1, flanking_is_zero)
  
  flanking_seg_is_cladewise_unreliable <- readRDS(file.path(folder, "processed", "flanking_seg_is_cladewise_unreliable", paste0(CHROM, ".rds")))
  flanking_pass <- (!flanking_seg_is_cladewise_unreliable)  & flanking_pass_2 #& no_unwanted_pass
  keep_after_first_t <- !assign_falses_after_the_first(!flanking_pass, !(flanking_is_zero > -1))
  is_enough <- sum(flanking_rep_len$group_count[keep_after_first_t]) >= min_supporting
  flanking_pass <- (flanking_pass * is_enough) == 1L
  
  is_cladewise_unreliable <- readRDS(file.path(folder, "processed", "is_cladewise_unreliable", paste0(CHROM, ".rds")))
  pass <- any(flanking_pass) & (!is_cladewise_unreliable) # & no_unwanted_l1 # & !l1_5p_is_unreliable
  
  # Tidy up the information
  reference_anno <- readRDS(file.path(folder, "processed", "reference_anno", paste0(CHROM, ".rds")))
  passed_anno <- reference_anno[pass]
  if(which_side == "upstream"){
    selected_info <- "fivep_frag"
    passed_flanking <- unlist(first_element((flanking_rep_len$rep_len[flanking_pass])[pass]))
    passed_anno$flanking_gcoord <- flank(resize(mcols(passed_anno)[[selected_info]], width = 1L, fix = "start"), width = passed_flanking) #to_check # 
    passed_anno$full_range <- unlist(range(split(c(empty_metadata(passed_anno), passed_anno$flanking_gcoord), seq_along(passed_anno))), use.names = FALSE)
    passed_anno$timepoint <- tp[pass]
  }
  
  out_folder <- paste0("predicted_chimera", "_", which_side)
  if(!dir.exists(file.path(folder, "processed", out_folder))){
    dir.create(file.path(folder, "processed", out_folder))
  }
  saveRDS(passed_anno, file.path(folder, "processed", out_folder, paste0(CHROM, ".rds")))
}
