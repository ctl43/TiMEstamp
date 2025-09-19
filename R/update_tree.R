#' Update Phylogenetic Tree
#'
#' @param tree A `treedata` object representing the phylogenetic tree.
#' @param species A vector of species to be kept, 
#'
#' @return A `treedata` object with updated metadata for sister clade groups.
#' @export
#' @importFrom ape keep.tip write.tree
update_tree <- function(tree, species, folder) {
  tree <- read.tree(tree)
  new_tree <- keep.tip(tree, species)
  write.tree(new_tree, file = file.path(folder, "updated_tree.nh"))
}
