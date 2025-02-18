#' Compute Alignment Coverage by Clades
#'
#' This function calculates the coverage of transposable element sequences across sister clades.
#' @param rmsk_range A GRanges object representing repeatMasker regions.
#' @param sister_clades A list of sister-clades.
#' @param folder A character string specifying the path to store output files.
#' @param folder A character string specifying the directory where output files will be stored.
#' @param chrom A character value specifying the chromosome of interest (default:NULL).
#' @export
get_coverage_by_clades <- function(rmsk_range, sister_clades, folder, chrom = NULL){
  if(is.null(chrom)){
    chrom <- sub(".rds$","",dir(file.path(folder, "full_gap")))
  }
  for(i in chrom){
    get_coverage_by_chrom(chrom = i, rmsk_range = rmsk_range, sister_clades = sister_clades, folder = folder)
  }
}

#' @importFrom IRanges Views coverage IntegerList ranges
#' @importFrom GenomicRanges GRanges setdiff
#' @importFrom GenomeInfoDb seqnames seqlevelsInUse seqlevels<-
get_coverage_by_chrom <- function(rmsk_range, sister_clades, chrom, folder){
  # Calculate missing coverage by clades
  missings <- readRDS(file.path(folder, "full_gap", paste0(chrom, ".rds")))
  missings_by_clades <- lapply(sister_clades, function(x){missings[missings$labels %in% x]})
  missing_cov_by_clades <- lapply(missings_by_clades, coverage)

  # Getting rmsk of a specific chr
  ref <- rmsk_range[seqnames(rmsk_range) == chrom]
  seqlevels(ref) <- seqlevelsInUse(ref)
  members_grp <- factor(rep(names(ref), lengths(ref$members)), levels = names(ref))

  # Extracting gap coverage
  extracted <- sapply(missing_cov_by_clades, function(x){
    val <- unlist(sapply(Views(x, unlist(ranges(ref$members))), sum), use.names = FALSE)
    sum(IntegerList(split(val, members_grp)))
  })

  # Getting total coverage * number of genomes in sister clades
  sum_size <- sapply(lengths(sister_clades), function(x){sum(x * width(ref$members))})
  missing_coverage <- extracted/sum_size
  aln_coverage <- 1 - missing_coverage
  if(!dir.exists(file.path(folder, "aln_coverage"))){
    dir.create(file.path(folder, "aln_coverage"))
  }
  saveRDS(aln_coverage, file.path(folder, "aln_coverage", paste0(chrom, ".rds")))
}
