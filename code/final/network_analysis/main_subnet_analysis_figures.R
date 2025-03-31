source("code/final/R_functions/taxa_names_functions.R")


source("code/final/network_analysis/get_spiec_adj_matrix.R")
source("code/final/network_analysis/get_spiec_edgelist.R")
source("code/final/network_analysis/filter_edglist_by_subnet_size.R")
source("code/final/network_analysis/create_subnet_figure_panel_A.R")

source("code/final/network_analysis/get_cor_hclust_order.R")

source("code/final/network_analysis/select_cox_phs_padj_for_vertices.R")

source("code/final/network_analysis/aaro_functions/create_cox_data.R")

source("code/final/network_analysis/run_cox_wrapper_with_covariates.R")
source("code/final/network_analysis/get_spiec_network_variables.R")

source("code/final/network_analysis/nice_otu_names_for_heat.R")
source("code/final/network_analysis/create_network_subheatmap.R")





  # run Cox on network item, create figures

main_subnet_analysis_figures <- function (pseq, datadir, pltdir, spiec.res, pseq_cox_res=NULL, fdr_levels=2) {


  text.size <-15
  legend.title.size <- 15
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


  ##### Cox regression (show Cox regression results as colors in network plots)
  # prelims for running survival analysis
  #pseq.cox <- pseq.rel.core.rel.clr

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


  # all functions to be made available in this repo.

  # obtain eschericiha component and plot it
  n.components <- subnet_components_from_edglist(filtd.edglist)
  n.components[["OTUsub"]] <- gsub(" \\(Bacteria\\)", "", n.components[["OTU"]])
  escherichia_component <- n.components[n.components[["OTUsub"]]=="Escherichia", "componentID"]
  vertices_in_escherichia_component <- filter(n.components, componentID==escherichia_component)[["OTU"]]
  escherichia.edglist <- filter(tmp.edglist,
    Var1 %in% vertices_in_escherichia_component | Var2 %in% vertices_in_escherichia_component)

  n_eschericiasn <- network::network(escherichia.edglist,
                                     directed=FALSE,
                                     matrix.type="edgelist",
                                     ignore.eval=FALSE)



  tmp <- set_up_n_and_palette(n_eschericiasn, pseq.core, cox_linear, fdr_levels=fdr_levels)
  n_eschericiasn <- tmp$n
  col_palette <- tmp$col_palette
  pA <- create_subnet_figure_panel_A(n_eschericiasn, col_palette)
  pA
  ggsave(paste0(pltdir, "snfig_esch_component_spieceasi_core_propcol_s", subnetthr*100,".png"), width=10, height=7)


  ##### Activation heatmap

  n.relabunds.core.genera.clrsc <- t(scale(t(abundances(pseq.rel.core.rel.clr))))

  n.relabunds.core.genera.clrsc.m <- melt(n.relabunds.core.genera.clrsc) %>% transmute(OTU=Var1,
                                                           SampleID=Var2,
                                                           scaled.CLR=value) %>% nice_otu_names_for_heat

  # order components by hclust

  transfd_abundances <- as.matrix(n.relabunds.core.genera.clrsc[vertices_in_escherichia_component,])
  order.components <- get_cor_hclust_order(t(transfd_abundances))

  # samples by relative abundances
  order.samples <- sample_names(pseq.rel.core)[order(colSums(abundances(pseq.rel.core)[vertices_in_escherichia_component,]))]

  tmp_fulln <- set_up_n_and_palette(n, pseq.core, cox_linear)
  n <- tmp_fulln$n

  pB <- create_network_subheatmap(filter(n.relabunds.core.genera.clrsc.m,
                                      OTU %in% vertices_in_escherichia_component),
                                  order.samples, order.components, n.components,
                                  value="scaled.CLR", n=n, text.size=text.size,
                                  legend.title.size=legend.title.size)
  pB

  ggsave(paste0(pltdir, "snfig_spieceasi_clusters_tmp", 100*subnetthr, ".png"), height=20, width=6)
  ggsave(paste0(pltdir, "snfig_spieceasi_clusters_tmp", 100*subnetthr, ".pdf"), height=20, width=6)


  ##### All together
  pall <- ggarrange(plot_grid(NULL, pB, rel_widths=c(0.04, 0.96), nrow=1),
                    plot_grid(NULL, pA, rel_widths=c(0.01, 0.99), nrow=1),
                    ncol=2, nrow=1, labels=c('A', 'B'), widths=c(1.5,1))

  pall
  ggsave(paste0(pltdir, "main.png"), width=16, height=4)
  ggsave(paste0(pltdir, "main.pdf"), width=16, height=4)

  return(list(pseq_cox_res, cox_linear, n.components)) # TODO add plot items pA, pB
}

# full figures, for supplements
