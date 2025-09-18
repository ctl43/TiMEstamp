#' Computing alignment coverage by clades
#'
#' This function calculates alignment coverage using repeat-masked (rmsk) regions
#' and phylogenetic tree information. It determines the sister-clades level alignment
#' coverage of the provided repeat annotation. Additionally, it will compute
#' if the region is a missing loci if specified.
#'
#' @param tree_file A character string specifying the path to the phylogenetic tree file.
#' @param rmsk_file A character string specifying the path to the repeat-masked annotation file.
#' @param folder A character string specifying the directory where output files will be stored.
#' @param consider_missing_loci A logical value indicating whether missing loci should be considered (default: TRUE).
#' @param chrom A character value specifying the chromosome of interest (default:NULL).
#'
#' @return No explicit return value; results are stored in the specified folder.
#' @export
get_alignment_coverage <- function(tree_file, rmsk_file, folder, consider_missing_loci = TRUE, chrom = NULL){
  if(is.null(rmsk_file)){
    rmsk_range <- readRDS(file.path(folder, "rmsk_range.rds"))
  }else{
    if(file.exists(rmsk_file)){
      rmsk_range <- import_rmsk(rmsk_file)
    }
  }
  sister_clades <- get_sister(tree_file)
  get_coverage_by_clades(rmsk_range = rmsk_range, sister_clades = sister_clades, folder = folder, chrom = chrom)
  if(consider_missing_loci){
    get_missing_loci(rmsk_range = rmsk_range, sister_clades = sister_clades, folder = folder, chrom = chrom)
  }
}
