
#' Identify missing loci by analyzing upstream and downstream gap lengths
#' if they exceed a specified threshold.
#'
#' This function determines missing loci by
#' computing the gap length for these missing loci using immediate and
#' buffered (±100bp) flanking regions to reduce potential alignment artifacts.
#'
#' The buffered gap length calculation prevents biases introduced by
#' artifacts at the very edge of the reference sequence by considering a
#' region 100bp upstream/downstream.
#'
#' @param rmsk_range A GRanges object representing repeatMasker regions.
#' @param sister_clades A list of sister-clades.
#' @param folder A character string specifying the path to store output files.
#' @param chrom A character value specifying the chromosome of interest (default:NULL).
#' @return Saves the computed gap lengths in `large_gap_len` and `buffered_large_gap_len directories` in the `folder` directory.
#' @export
get_missing_loci <- function(rmsk_range, sister_clades, folder, chrom = NULL){
  if(is.null(chrom)){
    chrom <- sub(".rds$","",dir(file.path(folder, "full_gap")))
  }
  for(i in chrom){
    get_missing_loci_by_chrom(chrom = i, rmsk_range = rmsk_range, sister_clades = sister_clades, folder = folder)
  }
}

#' @importFrom GenomeInfoDb seqnames seqlevelsInUse seqlevels<-
get_missing_loci_by_chrom <- function(rmsk_range, sister_clades, chrom, folder){
  ref <- rmsk_range[seqnames(rmsk_range) == chrom]
  seqlevels(ref) <- seqlevelsInUse(ref)
  missings <- readRDS(file.path(folder, "large_gap", paste0(chrom, ".rds")))
  missings <- missings[missings$labels %in% unlist(sister_clades)]
  missings_by_clades <- lapply(sister_clades, function(x){
    z <- missings[missings$labels %in% x]
    z$labels <- factor(as.character(z$labels), levels = x)
    z
  })
  if(!dir.exists(file.path(folder, "large_gap_len"))){
    dir.create(file.path(folder, "large_gap_len"))
  }
  gap_len <- lapply(missings_by_clades, compute_flanking_gap_length, ref = ref)
  saveRDS(gap_len, file.path(folder, "large_gap_len", paste0(chrom, ".rds")))

  if(!dir.exists(file.path(folder, "buffered_large_gap_len"))){
    dir.create(file.path(folder, "buffered_large_gap_len"))
  }
  buffered_gap_len_100 <- lapply(missings_by_clades, compute_flanking_gap_length, ref = ref + 100)
  saveRDS(buffered_gap_len_100, file.path(folder, "buffered_large_gap_len", paste0(chrom, ".rds")))
  rm(gap_len);rm(buffered_gap_len_100);gc()
}
