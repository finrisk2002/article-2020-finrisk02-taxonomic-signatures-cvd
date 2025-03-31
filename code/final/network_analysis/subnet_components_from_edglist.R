#' Return subnet components given an edge list
#'
#' "Subnet components" here means the components of an undirected network
#' with no edges between them.
#'
#' @param edglist a list of edges
#'
#' @return data.frame with columns componentID, OTU, OTUunderscored
#'
#' @export
#' ## NB importFrom does not work if this file is sourced
#' @importFrom igraph graph_from_edgelist 
#' @importFrom igraph components
subnet_components_from_edglist <- function(edglist) {

  n.ig <- igraph::graph_from_edgelist(as.matrix(edglist[,1:2]), directed=FALSE)

  n.components <- data.frame(componentID=igraph::components(n.ig)$membership)
  n.components[,"OTU"] <- rownames(n.components)
  n.components[,"OTUunderscored"] <- sapply(n.components[,"OTU"], taxa2underscore)
  n.components
}
