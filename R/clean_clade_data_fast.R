#' Clean clade-level missing coverage (fast path)
#'
#' Validates inputs produced by the fast clade pipeline and, for each requested
#' chromosome, filters gap intervals, computes (buffered) flanking-gap lengths,
#' identifies loci whose entire segments are absent, and masks clade-level
#' missing-coverage matrices where available data are insufficient.
#'
#' @param folder Character. Project folder containing \file{sister_clades.rds},
#'   \file{gap/}, and fast outputs under \file{fast/}.
#' @param min_gap_len Integer (bp). Discard gap intervals shorter than this
#'   before computing flanking-gap lengths (default \code{500}).
#' @param min_gap_len_for_abs Integer (bp). Threshold for calling a locus
#'   absent using flanking gaps (default \code{1000}).
#' @param buffered_size Integer vector (bp). Additional symmetric buffers used
#'   to recompute buffered flanking-gap lengths (default \code{c(100, 250)}).
#' @param min_available_data Numeric in \code{[0,1]}. Minimum fraction of
#'   species per clade required to keep a clade-level value; entries with fewer
#'   contributing species are set to \code{NA} (default \code{0.1}).
#' @param min_clade_size Integer. Clades with size \code{\le} this threshold are
#'   exempt from \code{min_available_data} masking (default \code{20}).
#' @param chrom Character vector or \code{NULL}. Chromosomes to process. If
#'   \code{NULL}, discovered from \file{folder/gap/*.rds}.
#'
#' @details
#' Expected fast inputs:
#' \itemize{
#'   \item \file{fast/reference_anno/CHROM.rds} — per-chromosome reference loci.
#'   \item \file{sister_clades.rds} — list of species per clade.
#'   \item \file{gap/CHROM.rds} — per-species gap \code{GRanges} with a \code{labels} factor.
#'   \item \file{fast/missings_by_clade/CHROM.rds} — clade-level missing coverage (as written by your pipeline).
#' }
#'
#' Outputs:
#' \itemize{
#'   \item \file{fast/cleaned_missing_by_clade/CHROM.rds} — clade-level matrix with
#'         low-coverage clade entries masked to \code{NA}.
#' }
#'
#' @return Invisibly returns \code{NULL}. Called for side effects (files written).
#'
#' @seealso \code{\link{clean_clade_data_by_loci_availability_by_chrom}}
#'
#' @export
clean_clade_data_fast <- function(folder, min_gap_len = 500, min_gap_len_for_abs = 1000, buffered_size = c(100, 250), min_available_data = 0.1, min_clade_size = 20, chrom = NULL){
  if(is.null(chrom)){
    chrom <- sub(".rds$","",dir(file.path(folder, "gap")))
  }
  
  if(!file.exists(file.path(folder, "sister_clades.rds"))){
    stop("No sister_clades file is found in the folder, please run get_missing_by_clades_fast first.")
  }
  
  reference_folder <- file.path(folder, "fast", "reference_anno")
  reference_file <- dir(reference_folder)
  if(file.exists(reference_folder)){
    if(!all(dir(reference_file) %in% paste0(chrom, ".rds"))){
      stop("No reference file is found in the folder, please run get_missing_by_clades_fast first.")
    }
  }
  
  if(is.null(chrom)){
    chrom <- sub(".rds$","",dir(file.path(folder, "gap")))
  }

  for(i in chrom){
    clean_clade_data_by_loci_availability_by_chrom(CHROM = i, folder = folder, min_gap_len = min_gap_len, min_gap_len_for_abs = min_gap_len_for_abs, 
                                                   buffered_size = buffered_size, min_available_data = min_available_data, min_clade_size = min_clade_size)
  }
}

#' @importFrom GenomeInfoDb seqnames seqlevelsInUse seqlevels<-
clean_clade_data_by_loci_availability_by_chrom <- function(CHROM, folder, min_gap_len, min_gap_len_for_abs, buffered_size, min_available_data, min_clade_size){
  # min_total_data: minimum number of genome assemblies required per clade to enable masking
  if(!dir.exists(folder)){
    stop("Output folder does not exist.")
  }
  reference_anno <- readRDS(file.path(folder, "fast", "reference_anno", paste0(CHROM, ".rds")))
  sister_clades <- readRDS(file.path(folder, "sister_clades.rds"))
  missings <- readRDS(file.path(folder, "gap", paste0(CHROM,".rds")))
  if(min_gap_len != 0){
    missings <- missings[width(missings) >= min_gap_len]; gc()
  }
  
  missings <- subsetByOverlaps(missings, range(reference_anno) + min_gap_len_for_abs + max(buffered_size) + 1000L)
  missings <- XVector::compact(missings)

  
  # Include the selected species only
  included_sp <- unlist(sister_clades, use.names = FALSE)
  missings <- missings[missings$labels %in% included_sp]
  missings$labels <- factor(missings$labels, levels = levels(missings$labels)[levels(missings$labels) %in% included_sp])
  
  # Computing the flanking gap of the L1
  flanking_gap_len <- list()
  flanking_gap_len[[1L]] <- compute_flanking_gap_length(missings, ref = reference_anno)
  
  if(length(buffered_size) != 0){
    for (i in seq_along(buffered_size)){
      flanking_gap_len[[i + 1L]] <- compute_flanking_gap_length(missings, ref = reference_anno + buffered_size[i])
    }
  }
  
  # Finding the species that miss the corresponding segment entirely
  is_abs_target_loci_final <- Reduce("|", lapply(flanking_gap_len, function(x)(x$up > min_gap_len_for_abs)|(x$dn > min_gap_len_for_abs)))
  
  # Clean up the data
  missing_by_clade <- readRDS(file.path(folder, "fast", "missing_by_clade", paste0(CHROM, ".rds")))
  
  # Assign NA if the available data per clade is less than 10% of the total assemblies
  if(min_available_data != 0){
    missing_by_clade <- mask_insufficient_clade_data(missing_by_clade, is_abs_target_loci_final, sister_clades, min_available_data = min_available_data, min_clade_size = min_clade_size)
  }
  if(!dir.exists(file.path(folder, "fast", "cleaned_missing_by_clade"))){
    dir.create(file.path(folder, "fast", "cleaned_missing_by_clade"))
  }
  saveRDS(missing_by_clade, file.path(folder, "fast", "cleaned_missing_by_clade", paste0(CHROM, ".rds")))
  
  rm(flanking_gap_len)
  rm(is_abs_target_loci_final)
  rm(missings);rm(reference_anno);gc()
}

mask_insufficient_clade_data <- function(missing_by_clade, loci_existence, sister_clades, min_available_data, min_clade_size){
  n_available_data_by_clade <- sapply(sister_clades, function(y)rowSums(!loci_existence[,y, drop = FALSE]))
  n_data_cutoff <- as.integer(lengths(sister_clades) * min_available_data)
  n_data_cutoff[lengths(sister_clades) <= min_clade_size]  <- 0 # If clades have less than 20 assemblies, do nothing
  n_data_cutoff <- t(replicate(nrow(n_available_data_by_clade), n_data_cutoff)) # can be improved
  missing_by_clade[n_available_data_by_clade < n_data_cutoff] <- NA
  return(missing_by_clade)
}
