#' Update and Save a Phylogenetic Tree
#'
#' Reads a phylogenetic tree, prunes it to retain only specified species, 
#' and writes the updated tree to disk in Newick format.
#'
#' @param tree Character or `phylo`. Input tree, either as a file path to a Newick-format tree 
#'   or as a `phylo` object (compatible with the `ape` package).
#' @param species Character vector. Names of species (tips) to retain in the tree.
#' @param folder Character. Output directory where the updated tree file 
#'   (`updated_tree.nh`) will be saved.
#'
#' @return A `phylo` object representing the pruned phylogenetic tree.
#'
#' @details
#' The function uses \code{\link[ape]{keep.tip}} to remove all taxa not listed in 
#' \code{species}, preserving the original branching structure among retained taxa.
#' The resulting tree is exported to \code{file.path(folder, "updated_tree.nh")}.
#'
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- rtree(10)
#' kept_species <- tree$tip.label[1:5]
#' pruned_tree <- update_tree(tree, kept_species, "results/")
#' }
#'
#' @export
#' @importFrom ape read.tree keep.tip write.tree
update_tree <- function(tree, species, folder) {
  tree <- read.tree(tree)
  new_tree <- keep.tip(tree, species)
  write.tree(new_tree, file = file.path(folder, "updated_tree.nh"))
  return(new_tree)
}
