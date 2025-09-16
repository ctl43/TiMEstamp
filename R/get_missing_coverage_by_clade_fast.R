#' Fast clade-level missing coverage
#'
#' Computes missing coverage per clade for selected chromosomes using
#' run-length coverage of per-species gap intervals. If \code{reference_file}
#' is not provided, a default \file{rmsk_range.rds} is used from \code{folder}.
#' Results are written under \file{fast/}.
#'
#' @param folder Character scalar. Project folder containing \file{gap/} and
#'   \file{sister_clades.rds}; outputs are written to \file{fast/}.
#' @param reference_file Character or \code{NULL}. Path to an RDS reference
#'   annotation. If \code{NULL}, defaults to \file{folder/rmsk_range.rds}.
#'   The object can be a \code{GRanges} or a \code{GRangesList} with a
#'   \code{$members} column (see Details).
#' @param chrom Character vector or \code{NULL}. Chromosomes to process.
#'   If \code{NULL}, discovered from \file{folder/gap/*.rds}.
#'
#' @details
#' \strong{Inputs expected in \code{folder}:}
#' \itemize{
#'   \item \file{sister_clades.rds}: list of character vectors of species per clade.
#'   \item \file{gap/CHROM.rds}: per-species gap \code{GRanges} with a factor
#'         column \code{labels} for species.
#'   \item \file{rmsk_range.rds} (optional): default reference if \code{reference_file} is \code{NULL}.
#' }
#'
#' \strong{Reference structure:}
#' \itemize{
#'   \item If the reference has a metadata column \code{$members} with a
#'         \code{GRangesList}, coverage is summed per member and then per locus.
#'   \item Otherwise a plain \code{GRanges} is used and coverage is summed per locus.
#' }
#'
#' \strong{Outputs written by the worker:}
#' \itemize{
#'   \item \file{fast/reference_anno/CHROM.rds}: per-chromosome subset of \code{reference_anno}.
#'   \item \file{fast/missing_by_clade/CHROM.rds}: numeric matrix (loci × clades)
#'         of missing coverage fractions in \[0,1].
#' }
#'
#' @return Invisibly returns \code{NULL}. Called for its side effects of creating
#'   the files listed above.
#'
#' @seealso \code{\link{get_missing_coverage_by_chrom}}
#'
#' @examples
#' \dontrun{
#' get_missing_coverage_by_clade_fast(folder = "proj")
#' get_missing_coverage_by_clade_fast(folder = "proj",
#'   reference_file = "proj/rmsk_range.rds", chrom = c("chr1","chr2"))
#' }
#'
#' @export
get_missing_coverage_by_clade_fast <- function(folder, reference_file = NULL, chrom = NULL){
  clade_file <- file.path(folder, "sister_clades.rds")
  if(!file.exists(clade_file)){
    stop("No clade information found. Please run get_sister() with your phylogenetic tree in Newick (NH) format first.")
  }
  sister_clades <- readRDS()
  
  if(is.null(reference_file)){
    reference_file <- file.path(folder, "rmsk_range.rds")
    if(!file.exists(reference_file)){
      stop("No default file is found in folder, please provide a valid reference file in RDS format.")
    }else{
      reference_anno <- readRDS(reference_file)
    }
  }
  
  if(is.null(chrom)){
    chrom <- sub(".rds$","",dir(file.path(folder, "gap")))
  }
  
  for(i in chrom){
    get_missing_coverage_by_chrom(chrom = i, reference_anno = reference_anno, sister_clades = sister_clades, folder = folder)
  }
}


#'
#' Internal worker that computes per-clade coverage over gaps for a given
#' chromosome, extracts coverage over the reference loci, and normalizes by
#' locus length to obtain missing coverage fractions.
#'
#' @param reference_anno GRanges or GRanges with a \code{$members} GRangesList column.
#' @param sister_clades List of character vectors; species per clade.
#' @param chrom Character scalar; chromosome identifier.
#' @param folder Project folder with \file{gap/} and \file{fast/}.
#'
#' @return Invisibly returns \code{NULL}. Writes:
#'   \file{fast/reference_anno/CHROM.rds} and
#'   \file{fast/missing_by_clade/CHROM.rds}.
#'
#' @keywords internal
#' @noRd
#' @importFrom GenomeInfoDb seqlevelsInUse seqinfo "seqinfo<-"
#' @importFrom GenomicRanges coverage GRanges ranges
#' @importFrom IRanges Views viewSums PartitioningByEnd IRanges width
get_missing_coverage_by_chrom <- function(reference_anno, sister_clades, chrom, folder){
  # Calculate missing coverage by clades
  missings <- readRDS(file.path(folder, "gap", paste0(chrom, ".rds")))
  seqinfo(missings) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)[chrom] # This has to be removed some day
  included_sp <- unlist(sister_clades, use.names = FALSE)
  missings <- missings[missings$labels %in% included_sp]
  missings$labels <- factor(missings$labels, levels = levels(missings$labels)[levels(missings$labels) %in% included_sp]);gc()
  missings <- split(missings, missings$labels)
  system.time(cov_by_species <- lapply(missings, function(x)coverage(x)[[1]]))
  rm(missings);gc()
  system.time(missing_cov_by_clades <- lapply(sister_clades, function(x){Reduce("+", cov_by_species[x])/length(x)}))
  rm(cov_by_species);gc()

  # Getting rmsk of a specific chr
  # Check what is inside the reference and hangle GRangeList and normal GRange differently
  has_members <- any(colnames(mcols(reference_anno)) == "members")
  reference_anno <- reference_anno[seqnames(reference_anno) == chrom]
  if(!dir.exists(file.path(folder, "fast", "reference_anno"))){
    dir.create(file.path(folder, "fast", "reference_anno"), recursive = TRUE)
  }
  if(has_members){
    ref <- reference_anno$members
    if(!methods::is(ref, "GRangesList")){
      stop("The reference_anno is not in expected structure.")
    }
    ref <- ref[lengths(ref) != 0L]
    seqlevels(ref) <- seqlevelsInUse(ref)
    sum_size <- sum(width(ref))
    members_grp <- PartitioningByEnd(ref)
    extracted <- sapply(missing_cov_by_clades, function(x){
      y <- Views(x, unlist(ranges(ref))) # ranges is not range!!!
      val <- viewSums(y)
      sum(split(val, members_grp))
    })
  }else{
    ref <- reference_anno
    if(!methods::is(ref, "GRanges")){
      stop("The reference_anno is not in expected structure.")
    }
    sum_size <- width(ref)
    extracted <- sapply(missing_cov_by_clades, function(x){
      y <- Views(x, ranges(ref)) # ranges is not range!!!
      viewSums(y)
    })
  }
  saveRDS(reference_anno, file.path(folder, "fast", "reference_anno", paste0(chrom, ".rds")))

  # Extracting gap coverage
  missing_coverage <- sweep(extracted, MARGIN = 1, STATS = sum_size, FUN = "/")
  if(!dir.exists(file.path(folder, "fast", "missing_by_clade"))){
    dir.create(file.path(folder, "fast", "missing_by_clade"))
  }
  saveRDS(missing_coverage, file.path(folder, "fast", "missing_by_clade", paste0(chrom, ".rds")))
}
