#' Estimates the insertion timepoint
#'
#' This function estimates the insertion timepoint of sequences based on
#' alignment coverage across sister clades, with an optional adjustment for missing loci.

#' @param folder A character string specifying the directory containing required data files.
#' @param consider_missing_loci A logical value indicating whether to account for missing loci (default: TRUE).
#' @param chrom A character string specifying the chromosome of interest.
#'
#' @return Saves the computed timepoint in `timepoint` in the specified `folder` directory.
#' @importFrom BiocGenerics match unlist
#' @export
get_timepoint <- function(folder, consider_missing_loci = TRUE, chrom = NULL){
  if(is.null(chrom)){
    chrom <- sub(".rds$","",dir(file.path(folder, "full_gap")))
  }

  for(i in chrom){
    get_timepoint_by_chrom(chrom = i, folder = folder, consider_missing_loci = consider_missing_loci)
  }
  rmsk_range <- readRDS(file.path(folder, "rmsk_range.rds"))
  original_id <- names(rmsk_range)
  rmsk_range <- rmsk_range[match(seqnames(rmsk_range), chrom, nomatch = 0) == 1L]
  rmsk_range <- clean_seqlevels(rmsk_range)
  rmsk_range <- split(rmsk_range, seqnames(rmsk_range))
  for(i in chrom){
    tp <- readRDS(file.path(folder, "timepoint", paste0(i, ".rds")))
    tp <- tp[names(rmsk_range[[i]])]
    rmsk_range[[i]]$timepoint <- tp
  }
  rmsk_range <- unlist(rmsk_range)
  rmsk_range$timepoint <- factor(rmsk_range$timepoint, levels = sort(unique(rmsk_range$timepoint)))
  saveRDS(rmsk_range, file.path(folder, "timepoint.rds"))
  gc()
}


#' Infer Insertion Timepoints by Chromosome
#' @param folder A character string specifying the directory containing required data files.
#' @param chrom A character string specifying the chromosome of interest.
#' @param consider_missing_loci A logical value indicating whether to account for missing loci (default: TRUE).
#'
#' @return Saves the computed timepoint in `timepoint` in the specified `folder` directory.
#' @export
get_timepoint_by_chrom <- function(folder, chrom, consider_missing_loci = TRUE) {
  # Load sister clades information
  sister_clades <- readRDS(file.path(folder, "sister_clades.rds"))

  # Load alignment coverage data for the specified chromosome
  aln_coverage <- readRDS(file.path(folder, "aln_coverage", paste0(chrom, ".rds")))

  if (consider_missing_loci) {
    # Load large gap length data
    gap_len <- readRDS(file.path(folder, "large_gap_len", paste0(chrom, ".rds")))

    # Load buffered large gap length data (±100bp extension)
    buffered_large_gap_len <- readRDS(file.path(folder, "buffered_large_gap_len", paste0(chrom, ".rds")))

    # Identify loci with large gaps in either dataset (upstream or downstream gaps)
    is_missing_locus <- mapply(
      function(x, y) x$up | x$dn | y$up | y$dn,
      x = gap_len, y = buffered_large_gap_len, SIMPLIFY = FALSE
    )
    is_missing_locus <- Map(function(x, y) x$up | x$dn | y$up | y$dn, gap_len, buffered_large_gap_len)

    # Count the number of loci that do not have large gaps
    n_with_loci <- sapply(is_missing_locus, function(x) rowSums(!x))

    # Determine data cutoff threshold (at least 10% of sister clades, min 5)
    n_data_cutoff <- as.integer(lengths(sister_clades) * 0.1)
    n_data_cutoff[n_data_cutoff <= 5] <- 0

    # Mask alignment coverage where has insufficient supporting data
    n_data_cutoff <- t(replicate(nrow(n_with_loci), n_data_cutoff))
    aln_coverage[n_with_loci < n_data_cutoff] <- NA
  }

  # Infer insertion timepoint using missing coverage (predict timepoint using missing coverage instead of alignment coverage)
  tp <- predict_timepoint(1 - aln_coverage)
  names(tp) <- rownames(aln_coverage)
  if(!dir.exists(file.path(folder, "timepoint"))){
    dir.create(file.path(folder, "timepoint"))
  }
  saveRDS(tp, file.path(folder, "timepoint", paste0(chrom, ".rds")))
}
