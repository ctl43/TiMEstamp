#' Estimate Insertion Timepoints from Clade-Level Coverage
#'
#' Estimates the evolutionary insertion timing of reference loci based on
#' clade-level coverage patterns. The function integrates per-chromosome
#' reference annotations with cleaned or raw missing-coverage matrices to infer
#' the earliest clade in which each locus is present.
#'
#' @param folder Character. Path to the project directory containing required
#'   subdirectories: \file{fast/reference_anno/}, \file{fast/missing_by_clade/}, 
#'   and optionally \file{fast/cleaned_missing_by_clade/}.
#' @param chrom Character vector or \code{NULL}. Chromosomes to process. If
#'   \code{NULL}, all chromosomes are inferred automatically from 
#'   \file{gap/*.rds}.
#'
#' @details
#' The function reads chromosome-specific reference annotations and corresponding
#' clade-level coverage matrices. If cleaned coverage data are available under
#' \file{fast/cleaned_missing_by_clade/}, they are used preferentially; otherwise,
#' the raw matrices in \file{fast/missing_by_clade/} are used with a warning.
#'
#' For each chromosome:
#' \enumerate{
#'   \item Load the reference loci from \file{fast/reference_anno/CHROM.rds}.
#'   \item Load the corresponding clade-level coverage matrix from
#'         \file{fast/cleaned_missing_by_clade/CHROM.rds} or
#'         \file{fast/missing_by_clade/CHROM.rds}.
#'   \item Predict the earliest clade of presence for each locus using
#'         \code{predict_timepoint()}.
#'   \item Append the predicted values to the reference object as a
#'         \code{timepoint} column.
#' }
#'
#' The merged results across chromosomes are saved as
#' \file{fast/predicted_tp.rds}, which contains a combined object (typically a
#' \code{GRanges}) with a \code{timepoint} factor column indicating the inferred
#' insertion clade index.
#'
#' @return Invisibly returns \code{NULL}. A combined \code{GRanges}-like object
#'   with predicted insertion timepoints is saved to
#'   \file{fast/predicted_tp.rds} within \code{folder}.
#'
#' @seealso
#' \code{\link{get_missing_by_clades_fast}},
#' \code{\link{clean_clade_data_fast}},
#' \code{\link{predict_timepoint}}
#'
#' @examples
#' \dontrun{
#' # Estimate insertion timepoints across all chromosomes
#' get_timepoint_fast(folder = "results/mammals_project")
#' }
#'
#' @importFrom BiocGenerics match unlist
#' @export

get_timepoint_fast <- function(folder, chrom = NULL){
  if(is.null(chrom)){
    chrom <- sub(".rds$","",dir(file.path(folder, "gap")))
  }

  reference_folder <- file.path(folder, "fast", "reference_anno")
  if(file.exists(reference_folder)){
    if(!all(paste0(chrom, ".rds") %in% dir(reference_folder))){
      stop("No reference file is found in the folder, please run get_missing_by_clades_fast on all chromosomes of interest first.")
    }
  }
  
  target_folder <- file.path(folder, "fast", "cleaned_missing_by_clade")
  if(dir.exists(target_folder)){
    if(!all(paste0(chrom, ".rds") %in% dir(target_folder))){
      stop("Not cleaned clade data is found in the folder, please run clean_clade_data_fast on all chromosomes of interest first.")
    }
  }else{
    warning("No cleaned clade data is found. Raw missing_by_clade was used.")
    target_folder <- file.path(folder, "fast", "missing_by_clade")
  }
  reference_files <- file.path(reference_folder, paste0(chrom, ".rds"))
  target_files <- file.path(target_folder, paste0(chrom, ".rds"))
  tp <- mapply(function(x, y){
    ref <- readRDS(x)
    target <- readRDS(y)
    ref$timepoint <- predict_timepoint(target)
    return(ref)
  }, x = reference_files, y = target_files, SIMPLIFY = FALSE)
  tp <- Reduce(c, tp)
  tp$timepoint <- factor(tp$timepoint, levels = min(tp$timepoint):max(tp$timepoint))
  saveRDS(tp, file.path(folder, "fast", "predicted_tp.rds"))
}
