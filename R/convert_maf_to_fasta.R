#' Convert all MAFs under a folder to FASTA, organized by chromosome
#'
#' @param maf_folder Character. Root folder containing MAF files (recursively).
#' @param species_file Character. Plain-text list of species IDs (one per line).
#' @param chrom_size_file Character. UCSC-style chrom.sizes for the reference assembly.
#' @param out_folder Character. Root output folder. Subfolders named by chromosome will be created.
#' @param reference_species Character. Optional. If empty, the reference is inferred from the first 's' line.
#' @param pattern Character. Regex to match MAF files. Default "^[^.].*\\.maf$" (plain .maf).
#' @param buffer_limit_mb Numeric. Memory buffer before flushing to disk.
#' @return Invisibly, the vector of processed MAF paths.
#' 
#' @importFrom BiocParallel bpmapply SnowParam
#' @export
convert_maf_to_fasta <- function(maf_folder,
                                 species_file,
                                 chrom_size_file,
                                 out_folder,
                                 reference_species = "",
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

  # Process each MAF; output goes to out_folder/fasta/<chr>/
  bpmapply(cxx_convert_maf_to_fasta,       
           maf_file          = maf_files,
           species_file      = species_file,
           chrom_size_file   = chrom_size_file,
           output_folder     = chrom_folders,
           reference_species = reference_species,
           buffer_limit_mb   = buffer_limit_mb,
           BPPARAM = SnowParam(workers = threads))

  invisible(maf_files)
}