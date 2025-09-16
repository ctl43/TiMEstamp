#' Extract gap runs from per-species FASTA files and write BED tracks per chromosome
#'
#' @description
#' Given a project `folder` that contains FASTA outputs produced by
#' `convert_maf_to_fasta()`, this function scans each chromosome subfolder under
#' `folder/fasta/`, detects contiguous gap runs (`-`) in each FASTA, and writes
#' one BED file per chromosome into `folder/gap_bed/`. The BED files contain
#' 0-based, half-open intervals of gap segments and include the FASTA filename
#' (species_tag) in the fourth column for provenance.
#'
#' @details
#' **Expected input layout**
#'
#' ```
#' <folder>/
#'   fasta/
#'     chr1/
#'       hg38.fa
#'       panTro6.fa
#'       gorGor6.fa
#'       ...
#'     chr2/
#'       ...
#'     chr22/
#'       ...
#' ```
#'
#' Each FASTA file contains a single record whose header is the reference
#' chromosome name (for example `>chr22`) and the sequence is the per-species
#' alignment string where unaligned positions are represented by `-`.
#'
#' **Outputs**
#'
#' For each `fasta/<chrom>/` subfolder, a BED file named `<chrom>.bed` is
#' written to `folder/gap_bed/`. Each BED row has four fields:
#'
#' 1. `chrom`: chromosome name taken from the FASTA header
#' 2. `start`: 0-based start coordinate of the gap run
#' 3. `end`: 0-based end coordinate (half-open)
#' 4. `name`: FASTA file basename without extension (species tag)
#'
#' Gap detection is performed by the compiled helper `cxx_fasta_gaps_to_bed()`,
#' which scans the sequences and emits gap intervals. By default, gaps of length
#' at least 10 bp are reported; adjust this threshold by passing
#' `min_gap_inclusive` to the C++ helper if needed.
#'
#' **Parallelization**
#'
#' Chromosome subfolders are processed in parallel with `BiocParallel::bpmapply`
#' using a `SnowParam` backend.
#'
#' @param folder Character scalar. Project root directory that contains a
#'   `fasta/` subdirectory with one subfolder per chromosome as produced by
#'   `convert_maf_to_fasta()`.
#' @param threads Integer scalar. Number of worker processes for parallel
#'   execution. Default is `1L` for serial execution.
#'
#' @return
#' Invisibly returns the vector of BED output file paths (one per chromosome).
#'
#' @section Requirements:
#' - The C++ routine `cxx_fasta_gaps_to_bed(fasta_files, bed_out, min_gap_inclusive)`
#'   must be available in the session (for example via `Rcpp::sourceCpp()`).
#' - FASTA files must contain the reference chromosome in the header and use
#'   `-` to indicate gaps.
#'
#' @seealso
#' `convert_maf_to_fasta()` for generating the per-species FASTA inputs;
#' `BiocParallel::bpmapply()`, `BiocParallel::SnowParam()`.
#'
#' @examples
#' \dontrun{
#' # Process gap tracks for all chromosomes under <folder>/fasta
#' extract_fasta_gaps_to_bed("path/to/project", threads = 4L)
#' }
#'
#' @importFrom BiocParallel bpmapply SnowParam
#' @importFrom rtracklayer import
#' @importFrom S4Vectors mcols
#' @export
extract_gaps_from_fasta <- function(folder, threads = 1L) {
  stopifnot(dir.exists(folder))
  fasta_folder <- file.path(folder, "fasta")   # produced by convert_maf_to_fasta
  stopifnot(dir.exists(fasta_folder))
  
  # chromosome subfolders under fasta/
  fa_folders <- dir(fasta_folder, full.names = TRUE)
  fa_folders <- fa_folders[file.info(fa_folders)$isdir]
  if (length(fa_folders) == 0L) stop("No chromosome subfolders found under: ", fasta_folder)
  
  # output folders
  gap_bed_folder <- file.path(folder, "gap_bed")
  if (!dir.exists(gap_bed_folder)) dir.create(gap_bed_folder, recursive = TRUE)
  
  gap_folder <- file.path(folder, "gap")
  if (!dir.exists(gap_folder)) dir.create(gap_folder, recursive = TRUE)

  # per-chromosome outputs
  bed_out <- file.path(gap_bed_folder, paste0(basename(fa_folders), ".bed"))
  gap_files <- file.path(gap_folder, paste0(basename(fa_folders), ".rds"))

  # list of FASTA files per chromosome (sorted)
  fasta_files <- lapply(fa_folders, function(x) {
    files <- dir(x, full.names = TRUE)
    files <- files[grepl("\\.fa(sta)?$", files, ignore.case = TRUE)]
    sort(files)
  })
  
  # Process each MAF; output goes to out_folder/<chr>/
  bpmapply(cxx_fasta_gaps_to_bed, 
           fasta_files = fasta_files,
           bed_out = bed_out,
           min_gap_inclusive = 10, 
           BPPARAM = SnowParam(workers = threads))
  
  gap_folder <- file.path(folder, "gap")
  if(!dir.exists(fasta_folder)){
    dir.create(gap_folder)
  }
  gap_files <- file.path(gap_folder, paste0(basename(fa_folders), ".rds"))
  mapply(function(bed, fa, out){
    gaps <- import(bed)
    sp <- sub(".fa$","",basename(fa))
    sp <- sort(sp)
    colnames(mcols(gaps)) <- "labels"
    mcols(gaps)[["labels"]] <- factor(mcols(gaps)[["labels"]], levels = sp)
    saveRDS(gaps, out)
  }, bed = bed_out, fa = fasta_files, out = gap_files)
}