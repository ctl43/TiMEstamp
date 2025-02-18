#' Update Phylogenetic Tree with New Sister Clades
#'
#' @param tree A `treedata` object representing the phylogenetic tree.
#' @param new_sister_clades A list where each element represents a sister clade, 
#'
#' @return A `treedata` object with updated metadata for sister clade groups.
#' @export
#' @importFrom tidytree full_join as_tibble as.treedata
update_tree <- function(tree, new_sister_clades) {
  # Create group for the sister clades
  list_grp <- as.character(rep(seq_along(new_sister_clades), lengths(new_sister_clades)))
  
  # Create a metadata dataframe linking sequence IDs to their corresponding groups
  meta_data <- data.frame(ID = unlist(new_sister_clades, use.names = FALSE), group = list_grp)
  
  # Merge metadata with the phylogenetic tree's existing data
  new_tree <- full_join(as_tibble(tree), meta_data, by = c("label" = "ID"))
  
  # Convert the updated tibble back into a `treedata` object
  new_tree <- as.treedata(new_tree)
  
  return(new_tree)
}
