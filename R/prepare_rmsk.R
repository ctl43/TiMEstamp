#' Import and Process RepeatMasker Output
#'
#' This function reads a RepeatMasker output file and processes it into a
#' `GRanges` object with metadata, organizing repeat elements by their IDs.
#'
#' @param rmsk_file A character string specifying the path to the RepeatMasker output file.
#' @param folder An optional character string specifying the directory to save the processed data. Default is NULL.
#'
#' @return A `GRanges` object containing processed repeat annotations with class, family, and repeat coordinates.
#' @export
#' @importFrom data.table fread
#' @importFrom S4Vectors DataFrame split mcols<-
#' @importFrom IRanges IRanges CharacterList
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocGenerics unlist unstrand
#' @importFrom XVector compact
prepare_rmsk <- function(rmsk_file, folder = NULL) {
  # Read RepeatMasker file, skipping first 3 lines and handling missing values
  rmsk <- fread(rmsk_file, fill = TRUE, skip = 3)
  rmsk <- DataFrame(rmsk)

  # Convert strand information: "C" to "-" for reverse strand
  rmsk[[9]][rmsk[[9]] == "C"] <- "-"

  # Adjust repeat start position based on strand orientation
  repeat_start <- rmsk[[12]]
  repeat_start[rmsk[[9]] == "-"] <- rmsk[[14]][rmsk[[9]] == "-"]
  repeat_start <- as.integer(repeat_start)

  # Adjust repeat left coordinate, removing parentheses
  repeat_left <- rmsk[[12]]
  repeat_left[rmsk[[9]] != "-"] <- rmsk[[14]][rmsk[[9]] != "-"]
  repeat_left <- as.integer(gsub("\\(|\\)", "", repeat_left))

  # Update processed start and left coordinates
  rmsk[[12]] <- repeat_start
  rmsk[[14]] <- repeat_left

  # Assign column names
  colnames(rmsk) <- c("score", "perc_div", "perc_del", "perc_ins", "chr", "start", "end", "left", "strand",
                      "repeat", "class/family", "repeat_start", "repeat_end", "repeat_left", "ID")

  # Remove invalid repeat elements where repeat_end is before repeat_start
  discarded_id <- rmsk[rmsk$repeat_end <= rmsk$repeat_start,]$ID
  rmsk <- rmsk[!rmsk$ID %in% discarded_id, ]

  # Parse into a `GRanges` object
  repeats <- GRanges(rmsk$chr, IRanges(rmsk$start, rmsk$end), strand = rmsk$strand)

  # Extract class and family from "class/family" column
  class_fam <- do.call(rbind, strsplit(rmsk$`class/family`, "\\/"))
  colnames(class_fam) <- c("class", "family")

  # Add metadata columns to the `GRanges` object
  mcols(repeats) <- cbind(rmsk[, c("repeat", "ID")], class_fam)

  # Store repeat coordinates separately
  repeats$repeat_coord <- GRanges(rmsk$`repeat`, IRanges(rmsk$repeat_start, rmsk$repeat_end))
  repeats$repeat_coord$left <- rmsk$repeat_left

  # Organize repeats by unique "repeat_ID" group
  rmsk_grp <- paste0(rmsk$`repeat`, "_", rmsk$ID)
  rmsk_grp <- factor(rmsk_grp, levels = unique(rmsk_grp))

  # Aggregate repeat sequences into CharacterLists
  rep_repeat <- CharacterList(split(rmsk$`repeat`, rmsk_grp))
  rmsk <- split(repeats, rmsk_grp)

  # Compute the range of each grouped repeat
  rmsk_range <- base::range(unstrand(rmsk))
  rmsk_range <- unlist(rmsk_range)

  # Assign metadata attributes to grouped repeats
  rmsk_range$`repeat` <- unlist(unique(rep_repeat), use.names = FALSE)
  rep_class <- CharacterList(split(repeats$class, rmsk_grp))
  rmsk_range$class <- unlist(unique(rep_class), use.names = FALSE)
  rep_family <- CharacterList(split(repeats$family, rmsk_grp))
  rmsk_range$family <- unlist(unique(rep_family), use.names = FALSE)
  rmsk_range$members <- rmsk

  # Remove empty elements
  rmsk_range <- compact(rmsk_range)

  # Optionally save the processed `GRanges` object
  if (!is.null(folder)) {
    saveRDS(rmsk_range, file.path(folder, "rmsk_range.rds"))
  }

  return(rmsk_range)
}
