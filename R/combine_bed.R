#' Combine and Process BED Files
#'
#' This function reads a metadata table containing BED file paths, chromosome information, 
#' and labels, then imports and processes the BED files into `GRanges` objects.
#' It stores the processed genomic ranges as `.rds` files for efficient access.
#' 
#' BED records exceeding a specified width (`large_gap_cutoff`) are stored separately.
#'
#' @param metatable A character string specifying the path to a metadata table (tab-separated format) 
#'   with columns: (1) BED file paths, (2) corresponding chromosomes, (3) labels (species or groups).
#' @param folder A character string specifying the directory where the output `.rds` files will be saved.
#' @param large_gap_cutoff An integer specifying the minimum width to classify as a "large gap" (default: 1000L).
#'
#' @return This function does not return an object but saves processed BED data as `.rds` files.
#'   - Full genomic ranges are saved in `"folder/full_gap/"`.
#'   - Large gaps are saved in `"folder/large_gap/"`.
#' @export
#' @importFrom data.table fread
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRangesList
combine_bed <- function(metatable, folder, large_gap_cutoff = 1000L) {
  # Read metadata table (tab-separated format)
  metatable <- data.table::fread(metatable)
  metatable <- as.data.frame(metatable)  # Convert to base R data frame for compatibility
  
  # Function to import BED files and process them into `GRanges`
  import_by_chrom <- function(x, sp, out, large_out) {
    # Import each BED file as a GRanges object
    grs <- GRangesList(lapply(x, rtracklayer::import))
    
    # Assign labels based on species or groups provided in `sp`
    labels <- factor(rep(sp, lengths(grs)), levels = sort(unique(sp)))
    
    # Flatten list into a single `GRanges` object
    gr <- unlist(grs, use.names = FALSE)
    gr$labels <- labels
    
    # Save full genomic ranges and large gaps separately
    saveRDS(gr, out)
    saveRDS(gr[width(gr) > large_gap_cutoff], large_out)
  }
  
  # Split metadata into chromosome-wise groups
  bed_by_chrom <- split(metatable[[1]], metatable[[2]])  # BED files grouped by chromosome
  sp_by_chrom <- split(metatable[[3]], metatable[[2]])   # Labels grouped by chromosome
  
  # Define output file paths
  out_by_chrom <- file.path(folder, "full_gap", paste0(names(bed_by_chrom), ".rds"))
  large_gap_out_by_chrom <- file.path(folder, "large_gap", paste0(names(sp_by_chrom), ".rds"))
  
  # Ensure output directories exist
  if (!dir.exists(file.path(folder, "full_gap"))) {
    dir.create(file.path(folder, "full_gap"))
  }
  if (!dir.exists(file.path(folder, "large_gap"))) {
    dir.create(file.path(folder, "large_gap"))
  }
  
  # Process each chromosome's data using `mapply`
  mapply(import_by_chrom, x = bed_by_chr, sp = sp_by_chrom, out = out_by_chrom, large_out = large_gap_out_by_chrom)
}
