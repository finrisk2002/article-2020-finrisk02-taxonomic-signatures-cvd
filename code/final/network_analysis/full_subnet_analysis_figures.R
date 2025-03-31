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



full_subnet_analysis_figures <- function (pseq, datadir, pltdir, spiec.res, pseq_cox_res=NULL, savefigs=TRUE,
                                          intermediate_pseqs=NULL, fdr_levels=2) {

  ### For supplements, we want a figure with the full network + activations


  text.size <-15
  # pretty names
  taxa_names(pseq) <- gsub("^g_", "", taxa_names(pseq))

  if (is.null(intermediate_pseqs)) {
    pseq.rel.core <- core(microbiome::transform(pseq, "compositional"),
      detection=0.1/100, prevalence=1/100)

    pseq.core <- prune_taxa(taxa_names(pseq.rel.core), pseq)
    pseq.rel.core.rel.clr <- transform(transform(pseq.rel.core, "compositional"), "clr")

    pseq.rel <- microbiome::transform(pseq, "compositional")
    pseq.rel.clr <- microbiome::transform(pseq.rel, "clr")
  } else {
    pseq.rel.core <- intermediate_pseqs$pseq.rel.core
    pseq.core <- intermediate_pseqs$pseq.core
    pseq.rel.core.rel.clr <- intermediate_pseqs$pseq.rel.core.rel.clr
    pseq.rel <- intermediate_pseqs$pseq.rel
    pseq.rel.clr <- intermediate_pseqs$pseq.rel.clr
  }

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
  

  tmp_fulln <- set_up_n_and_palette(n, pseq.core, cox_linear, fdr_levels=fdr_levels)
  n <- tmp_fulln$n
  col_palette <- tmp_fulln$col_palette

  pA_allsns <- create_subnet_figure_panel_A(n, col_palette)

  if (savefigs) {
    pA_allsns
    ggsave(paste0(pltdir, "snfig_spieceasi_core_propcol_s", subnetthr*100,".png"), width=10, height=7)
  }


  ### B panel
  ### Activation heatmap

  n.relabunds.core.genera.clrsc <- t(scale(t(abundances(pseq.rel.core.rel.clr))))

  n.relabunds.core.genera.clrsc.m <- melt(n.relabunds.core.genera.clrsc) %>% transmute(OTU=Var1,
                                                           SampleID=Var2,
                                                           scaled.CLR=value) %>% nice_otu_names_for_heat

  # it makes no sense to order this time by total abundances, so use hclust

  mat <- as.matrix(n.relabunds.core.genera.clrsc[rownames(n.relabunds.core.genera.clrsc) %in% (n %v% "taxa.names"),])
  order.samples <- get_cor_hclust_order(mat)

  n.components <- subnet_components_from_edglist(filtd.edglist)
  n.components[["OTUsub"]] <- gsub(" \\(Bacteria\\)", "", n.components[["OTU"]])
  # create a list of subheatmap plots for each subnet in n.components
  plist <- lapply(unique(n.components$componentID), function (cid) {
                                  genera.in.cid <- filter(n.components, componentID==cid)$OTU
                                  mat.sub <- as.matrix(n.relabunds.core.genera.clrsc[genera.in.cid,])
                                  # TODO: consider other beta div distance metrics here?
                                  order.components <- get_cor_hclust_order(t(mat.sub))
                                  create_network_subheatmap(filter(n.relabunds.core.genera.clrsc.m,
                                                                    OTU %in% genera.in.cid),
                                                             order.samples, order.components,
                                                             n.components,
                                                             value="scaled.CLR", n=n)
                                  })



  # heights parameters for plotting
  tmpheights <- sapply(unique(n.components$componentID), function(cid) {
                                                      1.8 + nrow(filter(n.components, componentID==cid))})
  pB_allsns <- ggarrange(plotlist=plist, nrow=length(plist), ncol=1, align="hv", 
                     legend="right", common.legend=TRUE,
                     labels=paste("SN", seq(length(plist))),
                     label.y=1.0, label.x=0.0, hjust=0.0, vjust=1.5,
                     heights=tmpheights, font.label=list(face="bold.italic", size=floor(0.85*text.size)))

  ### combine
  pall_allsns <- ggarrange(plot_grid(NULL, pB_allsns, rel_widths=c(0.04, 0.96), nrow=1),
                    plot_grid(NULL, pA_allsns, rel_widths=c(0.01, 0.99), nrow=1),
                    ncol=1, nrow=2, labels=c('A', 'B'), heights=c(2,2))

  if (savefigs) {
    pall_allsns
    ggsave(paste0(pltdir, "snfig_all_sns_acts_suppl.png"), width=10, height=15)
    ggsave(paste0(pltdir, "snfig_all_sns_acts_suppl.pdf"), width=10, height=15)
  }

  return(list(pseq_cox_res=pseq_cox_res,
            cox_linear=cox_linear,
            n.components=n.components,
            pA_allsns = pA_allsns,
            pB_allsns = pB_allsns))
}