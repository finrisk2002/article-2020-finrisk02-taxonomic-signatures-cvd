select_cox_phs_padj_for_vertices <- function(cox_linear, vertex_names) {

  tmp <- cox_linear %>% filter(predictor %in% vertex_names) %>% select(c("predictor","PH","p_adj"))
  tmp$vertex_cox_padj <- tmp$p_adj
  tmp$vertex_cox_phs <- tmp$PH

  return(tmp)
}
