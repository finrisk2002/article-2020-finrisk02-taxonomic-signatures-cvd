#' Return an adjacency matrix from a SpiecEasi computation result
#'
#' Note that lower trian. part of the matrix is set to 0.
#'
#' @param spiec.res SpiecEasi computation result
#' @param pseq  corresponding pseq object
#'
#' @return adjacency matrix
#'
#' @export
#'
#' @importFrom SpiecEasi getRefit
#' @importFrom microbiome abundances
get_spiec_adj_matrix <- function (spiec.res, pseq) {

  spiec.res.a <- as.matrix(getRefit(spiec.res))

  tmpnames <- my_pretty_network_vertex_names(pseq)
  colnames(spiec.res.a) <- rownames(spiec.res.a) <- colnames(t(abundances(pseq)))

  diag(spiec.res.a) <- 0
  spiec.res.a[lower.tri(spiec.res.a)] <- 0

  spiec.res.a
}