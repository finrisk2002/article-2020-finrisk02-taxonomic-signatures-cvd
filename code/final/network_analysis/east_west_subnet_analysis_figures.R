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



east_west_subnet_analysis_figure <- function (pseq,
 datadir, pltdir, spiec.res.east, spiec.res.west, pseq_cox_res=NULL, fdr_levels=2) {

  pseq.rel.core <- core(microbiome::transform(pseq, "compositional"),
      detection=0.1/100, prevalence=1/100)

  pseq.core <- prune_taxa(taxa_names(pseq.rel.core), pseq)
  pseq.rel.core.rel.clr <- transform(transform(pseq.rel.core, "compositional"), "clr")

  pseq.rel <- microbiome::transform(pseq, "compositional")
  pseq.rel.clr <- microbiome::transform(pseq.rel, "clr")


  # east west separation for all, not on core microbiome. done to retain all OTU names for refenes
  pseq.east <- prune_samples(sample_names(pseq)[meta(pseq)[,"EAST"] == "EAST"],
                                                 pseq)
  pseq.west <- prune_samples(sample_names(pseq)[meta(pseq)[,"EAST"] == "WEST"],
                                                 pseq)

  # helper
  s_east <- sample_names(pseq.east)
  s_west <- sample_names(pseq.west)

  # east west separtion on the core microbiome, the pseq analysis are done on.
  # not all intermediate_pseqs are used, but retained for coherence

  res_east <- full_subnet_analysis_figures(pseq.east, datadir, pltdir, spiec.res.east, savefigs=FALSE,
                  intermediate_pseqs=list(pseq.rel.core=prune_samples(s_east, pseq.rel.core),
                                          pseq.core=prune_samples(s_east, pseq.core),
                                          pseq.rel.core.rel.clr=prune_samples(s_east, pseq.core),
                                          pseq.rel=prune_samples(s_east, pseq.rel),
                                          pseq.rel.clr=prune_samples(s_east, pseq.rel.clr)),
                  fdr_levels=fdr_levels)
  
  res_west <- full_subnet_analysis_figures(pseq.west, datadir, pltdir, spiec.res.west, savefigs=FALSE,
                  intermediate_pseqs=list(pseq.rel.core=prune_samples(s_west, pseq.rel.core),
                                          pseq.core=prune_samples(s_west, pseq.core),
                                          pseq.rel.core.rel.clr=prune_samples(s_west, pseq.core),
                                          pseq.rel=prune_samples(s_west, pseq.rel),
                                          pseq.rel.clr=prune_samples(s_west, pseq.rel.clr)),
                  fdr_levels=fdr_levels)

  p_east <- res_east$pA_allsns
  p_west <- res_west$pA_allsns

  pew <- ggarrange(p_east, p_west, nrow=1, ncol=2,labels=c('EAST', 'WEST'))
  pew
  ggsave(paste0(pltdir, "snfig_east_west_suppl.png"), width=16, height=7)

  return(list(res_east=res_east, res_west=res_west, pew=pew))
}
