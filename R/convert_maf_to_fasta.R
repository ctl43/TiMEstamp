#' Convert all MAFs under a folder to FASTA, organized by chromosome
#'
#' @param maf_folder Character. Root folder that contains MAF files (searched recursively).
#' @param species_file Character or `NULL`. Path to a plain-text file listing
#'   species to retain (one per line). If `NULL`, species are inferred from
#'   `tree` and the list is written to `folder/species_list.txt`.
#' @param tree Character. Path to a Newick tree file used to curate the
#'   species list when `species_file` is `NULL`. Ignored when `species_file`
#'   is provided.
#' @param chrom_size_file Character. UCSC-style chrom.sizes for the reference assembly.
#' @param folder Character. Root output folder. Subfolders named by chromosome will be created.
#' @param pattern Character. Regular expression used to find MAF files.
#'   Default is `^[^.].*\\.maf$`, which matches plain `.maf` files and skips hidden files.
#' @param buffer_limit_mb Numeric. Memory buffer before flushing to disk.
#' @return Invisibly, the vector of processed MAF paths.
#' 
#' @importFrom BiocParallel bpmapply SnowParam
#' @importFrom ape read.tree
#' @export
convert_maf_to_fasta <- function(maf_folder,
                                 species_file = NULL, 
                                 tree = NULL,
                                 chrom_size_file,
                                 folder,
                                 pattern = "^[^.].*\\.maf$",
                                 buffer_limit_mb = 1L,
                                 threads = 1L) {
  # Getting full path
  maf_folder <- normalizePath(maf_folder, mustWork = TRUE)
  chrom_size_file <- normalizePath(chrom_size_file, mustWork = TRUE)
  folder <- normalizePath(folder, mustWork = FALSE)
  
  if (!dir.exists(folder)) dir.create(file.path(folder, "fasta"), recursive = TRUE, showWarnings = FALSE)
  
  if (!is.null(species_file)) {
    species_file <- normalizePath(species_file, mustWork = TRUE)
    
  } else if (!is.null(tree)) {
    tree <- normalizePath(tree, mustWork = TRUE)
    
    # Derive species file from the provided tree
    if (!file.exists(tree)) {
      stop("The specified tree file does not exist: ", tree)
    }
    
    tree_obj <- read.tree(tree)
    species_file <- file.path(folder, "species_list.txt")
    writeLines(tree_obj$tip.label, species_file)
    message("No species file provided. Generated a species list from the tree and saved to: ", species_file)
    
  } else {
    stop(
      "Either 'tree' or 'species_file' must be provided.\n",
      "Please specify a phylogenetic tree file (Newick format) or a text file containing species names."
    )
  }
  
  # Find MAFs recursively
  maf_files <- list.files(maf_folder, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(maf_files) == 0L) stop("No MAF files found under: ", maf_folder)
  chroms <- sub(".maf$", "", basename(maf_files))
  chrom_folders <- file.path(folder, "fasta", chroms)
  
  # Process each MAF; output goes to out_folder/fasta/<chr>/
  bpmapply(cxx_convert_maf_to_fasta,       
           maf_file          = maf_files,
           species_file      = species_file,
           chrom_size_file   = chrom_size_file,
           output_folder     = chrom_folders,
           reference_species = "",
           buffer_limit_mb   = buffer_limit_mb,
           BPPARAM = SnowParam(workers = threads))

  invisible(maf_files)
}