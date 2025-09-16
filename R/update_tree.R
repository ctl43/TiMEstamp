#' Update Phylogenetic Tree
#'
#' @param tree A `treedata` object representing the phylogenetic tree.
#' @param species A vector of species to be kept, 
#'
#' @return A `treedata` object with updated metadata for sister clade groups.
#' @export
#' @importFrom tidytree full_join as_tibble as.treedata
update_tree <- function(tree, species, folder) {
  tree <- read.tree(tree)
  new_tree <- keep.tip(tree, species)
  
  # # Create group for the sister clades
  # list_grp <- as.character(rep(seq_along(new_sister_clades), lengths(new_sister_clades)))
  # 
  # # Create a metadata dataframe linking sequence IDs to their corresponding groups
  # meta_data <- data.frame(ID = unlist(new_sister_clades, use.names = FALSE), group = list_grp)
  # 
  # # Merge metadata with the phylogenetic tree's existing data
  # new_tree <- full_join(as_tibble(tree), meta_data, by = c("label" = "ID"))
  # 
  # # Convert the updated tibble back into a `treedata` object
  # new_tree <- as.treedata(new_tree)
  write.tree(new_tree, file = file.path(folder, "updated_tree.nh"))
  # saveRDS(new_tree, file.path(folder, "sister_clades.rds"))
  # return(new_tree)
}
