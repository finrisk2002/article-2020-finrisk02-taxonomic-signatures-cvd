source("code/final/R_functions/taxa_names_functions.R")
source("code/final/network_analysis/set_ggnet_fdr2_col_pallette.R")
source("code/final/network_analysis/set_network_attributes_for_plots.R")
source("code/final/network_analysis/ggnet2repel.R")

# TODO move to create_...R
# create nice proportinal hazard ratio labels for ggnet2 plots
create_visual_ph_labels <- function (vertex_cox_phs,
                                     vertex_cox_padj,
                                     predictor,
                                     fdrthr = 0.05,
                                     fdrthr2 = 0.25) {
  cox_padj <- vertex_cox_padj

  cox.cuts <- cut(vertex_cox_phs, c(min(vertex_cox_phs), 1.0, max(vertex_cox_phs)),
                                        include.lowest=TRUE)

  # copy
  cox.cuts2 <-  cox.cuts

  not_sign_text <- paste("not signif.")
  not_sign_text2 <- paste("not signif.")

  levels(cox.cuts) <- c(levels(cox.cuts), not_sign_text)
  cox.cuts[cox_padj >= fdrthr] <- not_sign_text
  
  lvl_st_thr <- paste("HR \u2264 1, P <", fdrthr)
  lvl_lt_thr <- paste("HR > 1, P <", fdrthr)
  lvl_st_thr2 <- paste("HR \u2264 1, P <", fdrthr2)
  lvl_lt_thr2 <- paste("HR > 1, P <", fdrthr2)

  # set levels
  levels(cox.cuts) <- c(lvl_st_thr,
                        lvl_lt_thr,
                        not_sign_text)


  levels(cox.cuts2) <- c(lvl_st_thr,
                         lvl_lt_thr,
                         lvl_st_thr2,
                         lvl_lt_thr2,
                         not_sign_text2)
  # match levels with contents

  cox.cuts2[cox_padj >= fdrthr2] <- not_sign_text2
  cox.cuts2[cox_padj <= fdrthr2 & cox_padj > fdrthr] <- ifelse(vertex_cox_phs[cox_padj <= fdrthr2 & cox_padj > fdrthr] < 1,
                                        lvl_st_thr2,
                                        lvl_lt_thr2)

  cox.cuts <- ordered(cox.cuts, levels=c("HR > 1, P < 0.05", "HR ≤ 1, P < 0.05", "not signif."))
  cox.cuts2 <- ordered(cox.cuts2, levels=c("HR > 1, P < 0.05", "HR > 1, P < 0.25","HR ≤ 1, P < 0.05", "HR ≤ 1, P < 0.25","not signif."))

  return(list(cox_labels_fdr=cox.cuts,
              cox_labels_fdr2=cox.cuts2,
              predictor=predictor))
}


set_up_n_and_palette <- function(n, pseq, cox_linear, fdr_levels=2) {
	vertex_names <- sapply(n %v% "vertex.names", taxa2underscore)

	vertex_cox_params <- select_cox_phs_padj_for_vertices(cox_linear,
														  vertex_names)
	vis_ph <- create_visual_ph_labels(vertex_cox_params$vertex_cox_phs,
									  vertex_cox_params$vertex_cox_padj,
									  vertex_cox_params$predictor)

  if (fdr_levels ==2) {
	  col_palette <- set_ggnet_fdr2_col_palette(vis_ph)

	  n <- set_network_attributes_for_plots(n, vis_ph, my_pretty_network_vertex_names(pseq))
  } else {
    col_palette <- set_ggnet_fdr_col_palette(vis_ph)
    n <- set_network_attributes_for_plots(n, vis_ph, my_pretty_network_vertex_names(pseq), fdr_levels=1)
  }

	list(n=n, col_palette=col_palette)
}


# subnet figure
# n: network object
# pseq : pseq that contains the names of vertices
# cox_linear: linear association results from cox_wrapper
# 	required columns: p_adj

create_subnet_figure_panel_A <- function(n, col_palette) {
  # set up plotting

	legend.title.size <- 17
	text.size <-15
	taxa.text.size <-7


  	# plot
	pA <- ggnet2repel(n, label = TRUE, label.size=floor(taxa.text.size), color="Prop.Hazards",
		edge.color="color", edge.size="SpiecEasiEdgeW",
		color.palette=col_palette, color.legend="Prop.Hazards",
		legend.size=text.size,
		geom_text_or_label="label",
		segment.size  = 0.6, force=14,
		segment.color = "black") + theme(legend.title = element_text(size=legend.title.size)) +
		 guides(color=guide_legend("Risk of Death"))


}


# no color
create_subnet_figure_panel_A_only_net <- function(n, text.size =15, taxa.text.size = 13/3,
												  legend.title.size =17,
												  edge.size="SpiecEasiEdgeW", size=9) {


  	# plot
	pA <- ggnet2repel(n, label = TRUE, label.size=floor(taxa.text.size),
		edge.color="color", edge.size=edge.size,
		legend.size=text.size,
		size=size,
		geom_text_or_label="label",
		segment.size  = 0.6, force=14,
		segment.color = "black") + theme(legend.title = element_text(size=legend.title.size))

}