#' Identify Sister Clades from a Phylogenetic Tree
#'
#' This function extracts sister clades from a given phylogenetic tree.
#' It assumes that the first tip label in the tree represents the reference genome.
#'
#' @param tree A `phylo` object representing the phylogenetic tree, or a character string 
#'   specifying the file path to a Newick-format tree.
#' @param folder An optional character string specifying a directory to save the sister clades as an RDS file. Default is NULL.
#'
#' @return A `CharacterList` object where each element represents the sister clades for a given node.
#'   If `folder` is provided, the output is also saved as `"sister_clades.rds"` in the specified folder.
#' @export
#' @importFrom S4Vectors rev
#' @importFrom IRanges CharacterList
#' @importFrom ape read.tree nodepath extract.clade
get_sister <- function(tree, folder = NULL) {
  # Read the phylogenetic tree if a file path is provided
  if (!inherits(tree, "phylo")) {
    if (is.character(tree) && file.exists(tree)) {
      phylo_tree <- read.tree(tree)
    } else {
      stop("Invalid tree input: must be a 'phylo' object or a valid file path.")
    }
  } else {
    phylo_tree <- tree
  }
  
  # Assume the first tip label represents the reference genome
  ref_sp <- phylo_tree$tip.label[1]
  
  # Get the total number of internal nodes
  nnode <- phylo_tree$Nnode
  
  # Compute paths from each tip to the root
  nodepath_by_sp <- lapply(1:(nnode + 1), function(x) {
    nodepath(phy = phylo_tree, from = x, to = nnode + 2)[-1]
  })
  names(nodepath_by_sp) <- phylo_tree$tip.label
  
  # List all internal nodes
  all_nodes <- seq_len(nnode) + (nnode + 1)
  
  # Extract clades along the path from the reference species to the root
  clades <- lapply(rev(nodepath(phylo_tree, 1, nnode + 2)[-1]), extract.clade, phy = phylo_tree)
  
  # Retrieve tip labels for each clade
  clade_labs <- CharacterList(lapply(clades, function(x) x$tip.label))
  clade_labs <- c(clade_labs, CharacterList(ref_sp)) # Include reference species
  clade_labs <- rev(clade_labs) # Reverse order to align with hierarchy
  
  # Identify sister clades by excluding previously processed clades
  counter_labs <- clade_labs[-length(clade_labs)]
  counter_labs <- c(CharacterList(character()), counter_labs)
  sister_clades <- mapply(function(x, y) x[!x %in% y], x = clade_labs, y = counter_labs, SIMPLIFY = FALSE)
  sister_clades <- CharacterList(sister_clades)
  
  # Save results if a folder is specified
  if (is.null(folder)) {
    return(sister_clades)
  } else {
    saveRDS(sister_clades, file.path(folder, "sister_clades.rds"))
  }
}