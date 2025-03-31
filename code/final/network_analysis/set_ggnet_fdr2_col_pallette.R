# set 2-level color palatter FDR-adjusted p-values for ggnet2
set_ggnet_fdr2_col_palette <- function(vis_ph) {
  col_palette <- c("red", "salmon1","blue","skyblue2",   "gray")
  names(col_palette) <- levels(vis_ph$cox_labels_fdr2)
  return(col_palette)
}

set_ggnet_fdr_col_palette <- function(vis_ph) {
  col_palette <- c("red", "blue", "gray")
  names(col_palette) <- levels(vis_ph$cox_labels_fdr)
  return(col_palette)
}