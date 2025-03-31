nice_otu_names_for_heat <- function (net_abud_matrix) {

  net_abud_matrix['OTUnice'] <- gsub("\\(Bacteria\\)", " ", net_abud_matrix$OTU)
  net_abud_matrix['OTUnice'] <- gsub("\\(BacteriaPlasmid\\)", "\\(Plasmid\\)", net_abud_matrix$OTUnice)
  net_abud_matrix

}