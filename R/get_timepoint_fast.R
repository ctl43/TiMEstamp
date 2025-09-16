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
    target_folder <- file.path(folder, "fast", "missings_by_clade")
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
