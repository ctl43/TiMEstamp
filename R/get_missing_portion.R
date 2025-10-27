#' Compute per-species missing fractions for target loci
#'
#' Iterates over selected chromosomes and computes, for each locus in the
#' reference annotation, the fraction of alignment that is missing for each
#' species. Results are written as per-chromosome matrices under
#' \file{processed/missing_portion/}. This wrapper validates inputs, restricts
#' the reference annotation to requested chromosomes, and dispatches work to
#' \code{get_missing_portion_by_chrom()}.
#'
#' @param chrom Character vector or \code{NULL}. Chromosome names to process.
#'   If \code{NULL}, chromosomes are discovered from filenames under
#'   \file{folder/gap} by stripping the \code{.rds} suffix.
#' @param reference_file Character scalar. Path to an RDS file containing the
#'   reference annotation as a \code{GRanges} or \code{CompressedGRangesList}.
#'   The object may carry metadata columns used by \code{selected_info}.
#' @param selected_info Character or \code{NULL}. Optional selector of a
#'   metadata column in the reference annotation that defines the genomic
#'   intervals to analyze. Must be one of \code{"members"}, \code{"fivep_frag"},
#'   \code{"threep_frag"}, or \code{NULL}. When \code{NULL}, the top-level
#'   \code{GRanges} or \code{GRangesList} is used as is.
#' @param folder Character scalar. Project folder that contains the gap files
#'   under \file{gap/} and will receive outputs under \file{processed/}.
#'
#' @details
#' \strong{Input layout expected under \code{folder}:}
#' \itemize{
#'   \item \file{gap/CHROM.rds}: a \code{GRanges} of gap intervals per species for chromosome \code{CHROM}. 
#'         The \code{labels} column must be a factor giving species identifiers.
#'   \item \file{sister_clades.rds}: a list of character vectors defining species per clade.
#' }
#'
#' \strong{Reference annotation in \code{reference_file}:}
#' \itemize{
#'   \item A \code{GRanges} or \code{CompressedGRangesList} of loci. 
#'   \item If \code{selected_info = "members"}, a \code{GRangesList} column
#'         \code{$members} is used. 
#'   \item If \code{selected_info = "fivep_frag"} or \code{"threep_frag"},
#'         the corresponding \code{GRanges} column is used.
#' }
#'
#' \strong{Outputs written by this wrapper and its worker:}
#' \itemize{
#'   \item \file{processed/reference_anno/CHROM.rds}: the per-chromosome subset of the reference annotation.
#'   \item \file{processed/missing_portion/CHROM.rds}: a numeric matrix with rows as loci and columns as species
#'         giving missing fractions in \[0, 1].
#' }
#' If a chromosome contains more than 50,000 reference entries, the function
#' stops with a message. Consider subsetting the annotation or using the
#' clade-level pipeline (\code{\link{get_missing_coverage_by_clade_fast}},
#' \code{\link{clean_clade_data_fast}}, \code{\link{get_timepoint_fast}})
#' to reduce runtime and memory usage.
#' 
#' @return Invisibly returns \code{NULL}. This function is called for its side effects
#' of creating per-chromosome missing-fraction matrices.
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item Prepare \file{gap/CHROM.rds} files and \file{sister_clades.rds}.
#'   \item Prepare \code{reference_file} as an RDS of \code{GRanges} or \code{GRangesList}.
#'   \item Call \code{get_missing_portion()} to populate \file{processed/missing_portion/}.
#'   \item Downstream functions such as \code{classify_existence()} and \code{inspect_loci()}
#'         can then be applied.
#' }
#'
#' @examples
#' \dontrun{
#' # Using all chromosomes discovered in folder/gap
#' get_missing_portion(
#'   chrom = NULL,
#'   reference_file = "reference_anno.rds",
#'   selected_info = NULL,
#'   folder = "my_project"
#' )
#'
#' # Using a metadata column that holds 5' fragments
#' get_missing_portion(
#'   chrom = c("chr1", "chr2"),
#'   reference_file = "reference_anno.rds",
#'   selected_info = "fivep_frag",
#'   folder = "my_project"
#' )
#' }
#'
#' @export
get_missing_portion <- function(chrom = NULL, reference_file, selected_info = NULL, folder){
  if(!(is.null(selected_info)|selected_info %in% c("members", "fivep_frag", "threep_frag"))){
    stop("selected_info can only be either, members, fivep_frag, threep_frag or NULL")
  }
  
  
  if(!grepl("\\.rds$|\\.RDS", basename(reference_file))){
    stop("The reference_file must be in RDS format.")
  }else{
    reference_anno <- readRDS(reference_file)
  }
  
  
  if(!dir.exists(folder)){
    stop("The output folder is not existed.")
  }
  
  if(!dir.exists(file.path(folder, "processed"))){
    dir.create(file.path(folder, "processed"))
  }
  
  if(is.null(chrom)){
    chrom <- sub(".rds$","", dir(file.path(folder, "gap"), pattern = ".rds$"))
  }
  
  clade_file <- file.path(folder, "sister_clades.rds")
  if(!file.exists(clade_file)){
    stop("No clade information found. Please run get_sister() with your phylogenetic tree in Newick (NH) format first.")
  }
  
  if(selected_info == "members"){ # GRangesList
    reference_anno <- reference_anno[all(seqnames(reference_anno$members) %in% chrom)]
    anno_seqnames <- unlist(seqnames(reference_anno$members))
  }
  
  if(selected_info == "fivep_frag"){ # GRanges
    reference_anno <- reference_anno[seqnames(reference_anno$fivep_frag) %in% chrom]
    anno_seqnames <- seqnames(reference_anno$fivep_frag)
  }
  
  if(selected_info == "threep_frag"){
    reference_anno <- reference_anno[seqnames(reference_anno$fivep_frag) %in% chrom]
    anno_seqnames <- seqnames(reference_anno$threep_frag)
  }
  
  if(is.null(selected_info)){
    message("No metadata selected; using the GRanges/GRangesList object provided in reference_file as input.")
    anno_seqnames <- seqnames(reference_anno)
    reference_anno <- reference_anno[anno_seqnames %in% chrom]
  }
  
  
  # Ensure annotation entries are restricted to the expected chromosome(s)
  anno_seqname_count <- table(anno_seqnames)
  
  
  if(any(anno_seqname_count > 50000)){
    too_many <- paste(names(anno_seqname_count)[anno_seqname_count > 50000], collapse = " ")
    stop(paste0(
      "The following chromosome(s) have more than 50,000 entries: ",
      too_many,
      ". Please subset your data before processing, or use the clade-level method (faster)."
    ))
  }
  
  for(i in chrom){
    get_missing_portion_by_chrom(CHROM = i, reference_anno = reference_anno, selected_info = selected_info, folder = folder)
  }
}


#' Compute missing fractions for one chromosome
#'
#' Internal worker that subsets annotation to \code{CHROM}, filters gaps to clade species,
#' and writes the per-species missing-fraction matrix.
#'
#' @param reference_anno GRanges or GRangesList; reference annotation.
#' @param selected_info Character or \code{NULL}; selector of annotation column.
#' @param CHROM Character; chromosome identifier.
#' @param folder Project folder with \code{gap/} and \code{processed/}.
#'
#' @family missing-fraction
#' @keywords internal
#' @noRd
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols
#' @importFrom XVector compact
get_missing_portion_by_chrom <- function(reference_anno, selected_info, CHROM, folder){
  if(is.null(names(reference_anno))|any(duplicated(names(reference_anno)))){
    warning("The 'reference_anno' object does not contain unique names. Genomic coordinates will be used as fallback names.")
    if(any(duplicated(as.character(reference_anno)))){
      stop("The genomic coordinates are not unique.")
    }
    names(reference_anno) <- as.character(reference_anno)
  }
  sister_clades <- readRDS(file.path(folder, "sister_clades.rds"))
  reference_anno <- reference_anno[seqnames(reference_anno) == CHROM]
  reference_anno <- clean_seqlevels(reference_anno)
  mcols(reference_anno)[[selected_info]] <- clean_seqlevels(mcols(reference_anno)[[selected_info]])
  reference_anno <- XVector::compact(reference_anno);gc()
  
  if(!dir.exists(file.path(folder, "processed","reference_anno"))){
    dir.create(file.path(folder, "processed", "reference_anno"))
  }
  saveRDS(reference_anno, file.path(folder, "processed", "reference_anno", paste0(CHROM, ".rds")))
  
  ######################
  if(!dir.exists(folder)){
    dir.create(folder)
  }
  if(!dir.exists(file.path(folder, "processed", "missing_portion"))){
    dir.create(file.path(folder, "processed", "missing_portion"), recursive = TRUE)
  }
  
  missings <- readRDS(file.path(folder, "gap", paste0(CHROM,".rds")))
  
  # Include the selected species only
  included_sp <- unlist(sister_clades, use.names = FALSE)
  missings <- missings[missings$labels %in% included_sp]
  missings$labels <- factor(missings$labels, levels = levels(missings$labels)[levels(missings$labels) %in% included_sp]);gc()
  
  # Select useful missing data only
  if(is.null(selected_info)){
    anno_names <- names(reference_anno)
    message("No metadata selected; using the GRanges/GRangesList object provided in reference_file as input.")
  }else{
    anno_names <- names(reference_anno)
    reference_anno <- mcols(reference_anno)[[selected_info]]
  }
  
  missings <- subsetByOverlaps(missings, range(reference_anno) + 1000L)
  missing_portion <- compute_missing_portion(missings, reference_anno, name = anno_names)
  saveRDS(missing_portion, file.path(folder, "processed", "missing_portion", paste0(CHROM, ".rds")))
}



