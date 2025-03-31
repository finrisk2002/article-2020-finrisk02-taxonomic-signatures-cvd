source("code/final/R_functions/taxa_names_functions.R")


source("code/final/network_analysis/get_spiec_adj_matrix.R")
source("code/final/network_analysis/get_spiec_edgelist.R")
source("code/final/network_analysis/filter_edglist_by_subnet_size.R")
source("code/final/network_analysis/create_subnet_figure_panel_A.R")
source("code/final/network_analysis/select_cox_phs_padj_for_vertices.R")


source("code/final/network_analysis/get_cor_hclust_order.R")

source("code/final/network_analysis/aaro_functions/create_cox_data.R")

source("code/final/network_analysis/run_cox_wrapper_with_covariates.R")
source("code/final/network_analysis/get_spiec_network_variables.R")
source("code/final/network_analysis/nice_otu_names_for_heat.R")
source("code/final/network_analysis/create_network_subheatmap.R")



nocolor_subnet_analysis_figure <- function (pseq, datadir, pltdir, spiec.res, pseq_cox_res=NULL) {

  ### For supplements, we want a figure with the full network + activations


  text.size <-15
  # pretty names
  taxa_names(pseq) <- gsub("^g_", "", taxa_names(pseq))

  pseq.rel.core <- core(microbiome::transform(pseq, "compositional"),
   detection=0.1/100, prevalence=1/100)

  pseq.core <- prune_taxa(taxa_names(pseq.rel.core), pseq)
  pseq.rel.core.rel.clr <- transform(transform(pseq.rel.core, "compositional"), "clr")

  pseq.rel <- microbiome::transform(pseq, "compositional")
  pseq.rel.clr <- microbiome::transform(pseq.rel, "clr")

  # threshold for including edge in subnets
  # minimum size of subnet (independent network component) to be included in the edgelist
  subnetthr <- 0
  subnet.minsize <- 3

  network_vars <- get_spiec_network_variables(spiec.res, pseq.core, subnetthr, subnet.minsize)

  # adjacency matrix
  spiec.res.a <- network_vars$spiec.adjacency
  # edgelist format
  tmp.edglist <- network_vars$original.edglist
  filtd.edglist <- network_vars$filtd.edglist

  ##### Full filtered network

  n <- network::network(filtd.edglist,
                        directed=FALSE,
                        matrix.type="edgelist",
                        ignore.eval=FALSE)

  ### cox results

  # run cox_wrapper if pseq_cox_res is not provided
  pseq.cox <- pseq.rel.clr


  if (!exists("pseq_cox_res") || is.null(pseq_cox_res)) {

    pseq_cox_res <- run_cox_wrapper_with_covariates(pseq.cox,
      covariates = c("BL_AGE", "BMI", "MEN", "CURR_SMOKE",
        "PREVAL_DIAB", "SYSTM", "BL_USE_RX_L",
        "BP_TREAT"),
      alpha_level =0.05)
    
  }

  # extract relevant things from pseq_cox_res$results
  # hackish. better would be update the cox_wrapper
  cox_linear <- pseq_cox_res$results %>% filter(association=="linear")
  cox_linear$p_adj <- cox_linear$p %>% p.adjust(method="BH")
  ### a panel OK
  

  tmp_fulln <- set_up_n_and_palette(n, pseq.core, cox_linear)
  n <- tmp_fulln$n

  pA_allsns <- create_subnet_figure_panel_A_only_net(n, text.size=40, taxa.text.size=60, legend.title.size=45,
                                                     edge.size=6, size=30) # no colors.
  pA_allsns

  ggsave(paste0(pltdir, "snfig_nocolors_20200122.png"), width=40, height=30)
  ggsave(paste0(pltdir, "snfig_nocolors_20200122.pdf"), width=40, height=30)

  return(list(pseq_cox_res, cox_linear, n, pA_allsns))
}