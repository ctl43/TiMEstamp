#' Combine Unaligned Regions by Chromosome
#'
#' This function aggregates unaligned (missing/gap) regions across species for a specific chromosome.
#' It assumes RDS files (created by `get_unaligned_regions()`) are stored in a folder with filenames following
#' a naming pattern like `chr1_species_unaligned.rds`. The function supports multi-threaded reading
#' and outputs a single `GRanges` object with chromosome labels.
#'
#' @param CHROM Character. The chromosome name (e.g., `"chr1"`).
#' @param ref_genome Character. Reference genome name used for labeling (default `"hg38"`).
#' @param input_folder Character. Folder containing a folder named `"unaligned"` that has per-species unaligned `.rds` files.
#' @param output_folder Character. Folder to write the combined output. Defaults to `input_folder`.
#' @param selected_clades Character vector. Optional subset of species or clades to include.
#' @param threads Integer. Number of threads to use for parallel file loading (default `1L`).
#'
#' @import GenomicRanges
#' @import S4Vectors
#' @import BiocParallel
#' @importFrom BiocParallel SnowParam bplapply
#' @importFrom S4Vectors Rle
#' @return Saves a combined `.rds` file in a subfolder `unaligned_by_chr/` under `output_folder`.
#' @export
#'
#' @examples
#' \dontrun{
#' combined_unaligned_by_chrom("chr1", input_folder = "out/", threads = 4)
#' }

get_unaligned_regions <- function(CHROM, 
                                        ref_genome = "hg38",
                                        folder,
                                        selected_clades = NULL, 
                                        threads = 1L) {
  extract_unaligned_regions_by_species
  if (!dir.exists(output_folder)) {
    stop("Output folder does not exist.")
  }
  
  # Find relevant .rds files
  missing_files <- dir(file.path(folder, "unaligned"), pattern = paste0(CHROM, "_.*rds$"), full.names = TRUE)
  total_sp <- sub(".*_", "", gsub("chr._|_unaligned|.rds", "", basename(missing_files)))
  names(missing_files) <- total_sp
  missing_files <- unlist(missing_files, use.names = FALSE)
  
  # Split files for parallel reading
  core_grp <- cut(seq_along(missing_files),
                  breaks = as.integer(seq(1, length(missing_files), len = threads)),
                  include.lowest = TRUE)
  missing_files_split <- unname(split(missing_files, core_grp))
  
  # Parallel load
  cur_missings <- bplapply(missing_files_split, function(x) lapply(x, readRDS), BPPARAM = SnowParam(workers = threads))
  cur_missings <- unlist(cur_missings, recursive = FALSE)
  cur_missings <- GRangesList(cur_missings)
  
  if (any(duplicated(total_sp))) {
    stop("Duplicate species detected in filenames.")
  }
  
  if (!ref_genome %in% total_sp) {
    total_sp <- c(ref_genome, total_sp)
  }
  
  # Label and flatten
  labels <- Rle(names(cur_missings), lengths(cur_missings))
  cur_missings <- unlist(cur_missings, use.names = FALSE)
  cur_missings$labels <- factor(labels, total_sp)
  
  # Filter if needed
  if (!is.null(selected_clades)) {
    cur_missings <- cur_missings[cur_missings$labels %in% selected_clades]
    new_levels <- levels(cur_missings$labels)[levels(cur_missings$labels) %in% selected_clades]
    cur_missings$labels <- factor(as.character(cur_missings$labels), levels = new_levels)
  }
  
  # Save output
  unaligned_folder <- file.path(folder, "gap")
  if (!dir.exists(unaligned_folder)) {
    dir.create(unaligned_folder, recursive = TRUE)
  }
  
  output_file <- file.path(unaligned_folder, paste0(CHROM, ".rds"))
  saveRDS(cur_missings, output_file)
}

#' Identify Unaligned Regions from MAF and FASTA Files
#'
#' This function reads a multiple alignment format (MAF) file and a corresponding
#' FASTA file to identify regions that are missing (`*`) or gaps (`-`) in the alignment,
#' relative to a reference genome. It returns merged genomic regions meeting a
#' minimum length requirement.
#'
#' @param fa Character. Path to the input FASTA file.
#' @param maf Character. Path to the input MAF file (only the first 3 lines are used).
#' @param ref_genome Character. Reference genome identifier (default is `"hg38"`).
#' @param out_folder Character. Output directory to save merged region results.
#' @param min_len Integer. Minimum length (in base pairs) of the merged region to be saved (default is `10L`).
#'
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @return Saves a merged `.rds` file of `GRanges` objects representing the gaps and missing regions.
#' @export
#'
#' @examples
#' \dontrun{
#' get_unaligned_regions("example.fa", "example.maf", "hg38", "output_folder")
#' }

extract_unaligned_regions_by_species <- function(fa, maf, ref_genome = "hg38", folder, min_len = 10L) {
  out <- sub(".fa$", "", basename(fa))
  fa <- readLines(fa)
  maf <- readLines(maf, n = 3)
  
  target <- grep(paste0("s ", ref_genome), maf, value = TRUE)
  target <- read.table(text = target, col.names = FALSE)[, 1]
  
  if (length(target) == 0) {
    stop("Reference genome not found in MAF file.")
  }
  
  chr <- sub(paste0(ref_genome, "."), "", target[2])
  offset <- as.integer(target[3])
  
  headers <- grep(">", fa)
  cur_seq <- paste(fa[(headers[2] + 1L):length(fa)], collapse = "") # The first header is the reference sequence, so it starts with the second sequence.
  
  missing <- gregexpr("[*]+", cur_seq)
  if (!is.na(missing[[1]][1]) && missing[[1]][1] != -1) {
    missing <- GRanges(chr, IRanges(start = missing[[1]] + offset, width = attr(missing[[1]], "match.length")))
  }
  
  gap <- gregexpr("[-]+", cur_seq)
  if (!is.na(gap[[1]][1]) && gap[[1]][1] != -1) {
    gap <- GRanges(chr, IRanges(start = gap[[1]] + offset, width = attr(gap[[1]], "match.length")))
  }
  
  unaligned_folder <- file.path(folder, "unaligned")
  if (!dir.exists(folder)) {
    stop("Output folder does not exist.")
  }
  if (!dir.exists(unaligned_folder)) {
    dir.create(unaligned_folder)
  }
  
  if (inherits(gap, "GRanges") && inherits(missing, "GRanges")) {
    combined <- merge_granges(c(gap, missing), tol = 0)$regions
    combined <- combined[width(combined) >= min_len]
    # Add a step to include seqinfo in the GRanges objects
    saveRDS(combined, file.path(unaligned_folder, paste0(out, "_unaligned.rds")), compress = TRUE)
  }
}
