source("code/final/network_analysis/subnet_components_from_edglist.R")

#' Remove edges belonging to small subnets
#'
#' Retrieve components (subnet_components_from_edglist), remove edges involving
#' components smaller than subnet.minsize.
#'
#' @param edglist list of edges
#' @param subnet.minsize minumum size for component
#'
#' @return thresholded edgelist
#'
#' @export
#'
#' @importFrom dplyr filter
filter_edglist_by_subnet_size <- function(edglist, subnet.minsize=0) {

  # filter by subnet size
  n.components <- subnet_components_from_edglist(edglist)

  # sum for each component
  cids <- unique(n.components$componentID)
  n.comp.sums <- sapply(cids, function (cid) { nrow(filter(n.components, componentID==cid)) })
  names(n.comp.sums) <- cids

  n.comp.sums.thrd <- n.comp.sums[n.comp.sums >= subnet.minsize]
  n.components.thrd <- filter(n.components, componentID %in% names(n.comp.sums.thrd))

  # filter
  edglist.thrd <- filter(edglist, Var1 %in% n.components.thrd$OTU | Var2 %in% n.components.thrd$OTU)

  edglist.thrd
}