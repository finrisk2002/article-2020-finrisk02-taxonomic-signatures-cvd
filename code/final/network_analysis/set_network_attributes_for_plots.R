# add ggnet2 plotting attributes to network::network object
set_network_attributes_for_plots <- function(n,
                                             vis_ph,
                                             vertex_plot_names,
                                             edge_color_name="SpiecEasiEdge",
                                             fdr_levels=2) {

  # edge color: sign of the correlation

  network::set.edge.attribute(n, "color",
                              ifelse(n %e% edge_color_name > 0, "darkred", "darkblue"))

  # vertex color: Cox model color, depending on FDR threshold
  # NOTE: continuous color scales are difficult with ggnet2, so we categorize the value range

  # TODO: prevalence as size


  # Cox PH HRs as colors, shade indicates FDR
  if (fdr_levels==2) {
    fdr_labels <- as.vector(as.character(vis_ph$cox_labels_fdr2))
  } else {
    fdr_labels  <- as.vector(as.character(vis_ph$cox_labels_fdr))
  }
  # ensure they are ordered correctly
  names(fdr_labels) <- vis_ph$predictor
  fdr_labels_ordered <- fdr_labels[sapply(n %v% "vertex.names", taxa2underscore)]
  network::set.vertex.attribute(n, "Prop.Hazards", fdr_labels_ordered)
  network::set.vertex.attribute(n, "taxa.names", (n %v% "vertex.names"))
  network::set.vertex.attribute(n, "vertex.names", 
                                gsub("_", " ",
                                     vertex_plot_names[n %v% "vertex.names"]))
  return(n)
}