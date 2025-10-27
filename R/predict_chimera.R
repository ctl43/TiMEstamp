#' Predict Chimeric Insertions
#'
#' Identifies loci with putative chimeric insertions by comparing the inferred
#' insertion timepoints of each locus and its flanking segments across sister
#' clades. Loci are called as chimeric if their flank and locus timepoints
#' differ while passing reliability, length, and support filters.
#'
#' @param folder Character scalar. Path to the project folder containing all
#'   required processed inputs and where outputs will be written under
#'   \file{processed/}.
#' @param chrom Character vector or \code{NULL}. Chromosomes to process. If
#'   \code{NULL}, all available chromosomes are inferred automatically from
#'   filenames in \file{folder/gap/}.
#' @param which_side Character string. Either \code{"upstream"} (5' flank) or
#'   \code{"downstream"} (3' flank) to specify which flanking region to evaluate.
#'   Default: \code{"upstream"}.
#' @param min_len Integer (bp). Minimum representative flank length to include
#'   as a candidate segment. Default: \code{25}.
#' @param max_len Integer (bp). Maximum representative flank length to include.
#'   Default: \code{1000}. The function aborts if \code{max_len < 1000} or if
#'   \code{max_len <= min_len}.
#' @param min_cov Numeric in \code{[0, 1]}. Minimum clade-level coverage agreement
#'   used internally as a missingness cutoff (\code{1 - min_cov}) for
#'   \code{predict_timepoint()}. Default: \code{0.65}.
#' @param min_supporting Integer. Minimum total number of supporting flanking
#'   segments (after filtering) required for a locus to be called. Default:
#'   \code{15L}.
#'
#' @details
#' The function orchestrates per-chromosome analysis via
#' \code{predict_chimera_by_chrom()}, performing the following steps:
#' \enumerate{
#'   \item Compute insertion timepoints for loci and flanking segments using
#'         \code{predict_timepoint()} with a coverage cutoff of \code{1 - min_cov}.
#'   \item Enforce segment length filters (\code{min_len}, \code{max_len}).
#'   \item Apply clade-level reliability filters to exclude unreliable loci and segments.
#'   \item Aggregate results across flanking segment clusters, requiring at least
#'         \code{min_supporting} total segment support per locus.
#'   \item Output a subset of reference loci meeting all chimera criteria.
#' }
#'
#' \strong{Expected Inputs (in \file{processed/}):}
#' \itemize{
#'   \item \file{cleaned_missing_by_clade/CHROM.rds} — clade-level missing fractions for loci.
#'   \item \file{cleaned_flanking_seg_missing_by_clade/CHROM.rds} — clade-level missing fractions for flanking segments.
#'   \item \file{flanking_rep_len/CHROM.rds} — representative flank lengths and clustering metadata.
#'   \item \file{flanking_seg_grp/CHROM.rds} — factor mapping segments to parent loci.
#'   \item \file{flanking_seg_is_cladewise_unreliable/CHROM.rds} — segment-level reliability flags.
#'   \item \file{is_cladewise_unreliable/CHROM.rds} — locus-level reliability flags.
#'   \item \file{reference_anno/CHROM.rds} — per-chromosome reference loci.
#' }
#'
#' \strong{Outputs:}
#' \itemize{
#'   \item \file{processed/predicted_chimera_<which_side>/CHROM.rds} — loci passing
#'         chimera criteria, with additional fields:
#'         \code{$flanking_gcoord} (flanking genomic coordinates) and
#'         \code{$full_range} (combined locus range).
#'   \item \file{processed/predicted_chimera_<which_side>.txt} — combined text
#'         summary of all chromosomes.
#' }
#'
#' Parameter constraints are validated internally: \code{min_len > 10},
#' \code{max_len >= 1000}, and \code{max_len > min_len}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effects (files written).
#'
#' @seealso
#' \code{\link{extract_flanking_segment}},
#' \code{\link{predict_timepoint}},
#' \code{\link{predict_chimera_by_chrom}}
#'
#' @examples
#' \dontrun{
#' # Predict upstream (5') chimeric insertions across all chromosomes
#' predict_chimera(folder = "results/mammals_project", which_side = "upstream")
#'
#' # Predict downstream (3') chimeric insertions for chromosome 1 only
#' predict_chimera(folder = "results/mammals_project", chrom = "chr1", which_side = "downstream")
#' }
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
