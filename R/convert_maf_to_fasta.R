#' Convert all MAFs under a folder to FASTA, organized by chromosome
#'
#' @param maf_folder Character. Root folder that contains MAF files (searched recursively).
#' @param species_file Character or `NULL`. Path to a plain-text file listing
#'   species to retain (one per line). If `NULL`, species are inferred from
#'   `tree_file` and the list is written to `out_folder/species_list.txt`.
#' @param tree_file Character. Path to a Newick tree file used to curate the
#'   species list when `species_file` is `NULL`. Ignored when `species_file`
#'   is provided.
#' @param chrom_size_file Character. UCSC-style chrom.sizes for the reference assembly.
#' @param out_folder Character. Root output folder. Subfolders named by chromosome will be created.
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
                                 out_folder,
                                 pattern = "^[^.].*\\.maf$",
                                 buffer_limit_mb = 1L,
                                 threads = 1L) {
  stopifnot(dir.exists(maf_folder))
  if (!dir.exists(out_folder)) dir.create(file.path(out_folder, "fasta"), recursive = TRUE, showWarnings = FALSE)
  
  # Find MAFs recursively
  maf_files <- list.files(maf_folder, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(maf_files) == 0L) stop("No MAF files found under: ", maf_folder)
  chroms <- sub(".maf$", "", basename(maf_files))
  chrom_folders <- file.path(out_folder, "fasta", chroms)
  
  if (is.null(tree) && is.null(species_file)) {
    stop("Either 'tree' or 'species_file' must be provided.\n",
         "Provide a phylogenetic tree file (Newick format) or a text file containing species names.")
  }
  
  # If no species file is provided, derive one from the tree
  if (is.null(species_file)) {
    if (!file.exists(tree)) {
      stop("The specified tree file does not exist: ", tree)
    }
    tree <- read.tree(tree)
    species_file <- file.path(out_folder, "species_list.txt")
    writeLines(tree$tip.label, species_file)
    message("No species file provided. Generated a species list from the tree and saved to: ", species_file)
  }
  
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