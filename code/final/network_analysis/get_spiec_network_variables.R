get_spiec_network_variables <- function(spiec.res, pseq.core, subnetthr, subnet.minsize) {

  # adjacency matrix
  spiec.res.a <- get_spiec_adj_matrix(spiec.res, pseq.core)

  # to edgelist format

  tmp.edglist <- get_spiec_edgelist(spiec.res, spiec.res.a, pseq.core, subnetthr,
      edgeWcoef=5)

  # deprecated function name
  #tmp.edglist <- my_spiec_edglist_with_trappings(spiec.res, spiec.res.a, pseq.core, subnetthr,
  #  edgeWcoef=5)
  filtd.edglist <- filter_edglist_by_subnet_size(tmp.edglist, subnet.minsize)

  list(spiec.adjacency=spiec.res.a, original.edglist=tmp.edglist,
       filtd.edglist=filtd.edglist)

}