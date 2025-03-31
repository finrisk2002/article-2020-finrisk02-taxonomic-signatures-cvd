# compare SpiecEasi results for antibiotics / without antibiotics

source("code/final/R_functions/isinstalled.R")
source("code/final/R_functions/taxa_names_functions.R")


source("code/final/network_analysis/get_spiec_adj_matrix.R")
source("code/final/network_analysis/get_spiec_edgelist.R")
source("code/final/network_analysis/get_spiec_network_variables.R")
source("code/final/network_analysis/filter_edglist_by_subnet_size.R")



is.installed(c("biomformat", "microbiome", "phyloseq", "data.table",
               "RColorBrewer", "dplyr", "reshape2", "SpiecEasi"))

source("code/final/network_analysis/pseq_spieceasi.R")

genus_level_phyloseq_path_woutab <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21_nop_norxj01.rds"

genus_level_phyloseq_path_withab <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21_nop.rds"

phfg.pseq.d50kwithab <- readRDS(genus_level_phyloseq_path_withab)

pseq.d50kwithab.rel.core <- core(microbiome::transform(phfg.pseq.d50kwithab, "compositional"),
                           detection=0.1/100, prevalence=1/100)

pseq.d50kwithab.core <- prune_taxa(taxa_names(pseq.d50kwithab.rel.core), phfg.pseq.d50kwithab)


phfg.pseq.d50kwoutab <- readRDS(genus_level_phyloseq_path_woutab)

pseq.d50kwoutab.rel.core <- core(microbiome::transform(phfg.pseq.d50kwoutab, "compositional"),
                           detection=0.1/100, prevalence=1/100)

pseq.d50kwoutab.core <- prune_taxa(taxa_names(pseq.d50kwoutab.rel.core), phfg.pseq.d50kwoutab)


datadir_withab <- "output_2020-04-06_wab/raw_output/"
datadir_woutab <- "output_2020-04-01/raw_output/"
spieceasi.res.th001_withab <- readRDS(paste0(datadir_withab,
    "/spieceasi_pseq_2018-12-21_drop50k_core_genus_glasso_nl30_th001.RDs"))

spieceasi.res.th001_woutab <- readRDS(paste0(datadir_woutab,
    "/spieceasi_pseq_2018-12-21_drop50k_core_genus_glasso_nl30_th001.RDs"))

subnetthr <- 0


network_vars_min0_withab <- get_spiec_network_variables(spieceasi.res.th001_withab, pseq.d50kwithab.core, subnetthr, subnet.minsize=0)


network_vars_min0_woutab <- get_spiec_network_variables(spieceasi.res.th001_woutab, pseq.d50kwoutab.core, subnetthr, subnet.minsize=0)

spiec.res.a0.withab <- network_vars_min0_withab$spiec.adjacency

spiec.res.a0.woutab <- network_vars_min0_woutab$spiec.adjacency


# do pseq.cores actually match?

ntaxa(pseq.d50kwoutab.core)
# [1] 91

ntaxa(pseq.d50kwithab.core)
# [1] 91

taxa_names(pseq.d50kwoutab.core)[!taxa_names(pseq.d50kwoutab.core) %in% taxa_names(pseq.d50kwithab.core)]
#[1] "Leuconostoc (Bacteria)"
taxa_names(pseq.d50kwoutab.core)[!taxa_names(pseq.d50kwithab.core) %in% taxa_names(pseq.d50kwoutab.core)]
#[1] "Anaerotruncus (Bacteria)"

taxa_common <- taxa_names(pseq.d50kwoutab.core)[taxa_names(pseq.d50kwoutab.core) %in% taxa_names(pseq.d50kwithab.core)]

spiec.res.a0.withab.cm <- spiec.res.a0.withab[taxa_common, taxa_common]
spiec.res.a0.woutab.cm <- spiec.res.a0.woutab[taxa_common, taxa_common]

# edges
sum(spiec.res.a0.withab.cm)
sum(spiec.res.a0.woutab.cm)

# edgelist
network_vars_min0_woutab$original.edglist[, c("Var1", "Var2")]
network_vars_min0_withab$original.edglist[, c("Var1", "Var2")]

# check
helper_edges  <- function (mymat) {
ni <- c()
nj <- c()
for (i in 1:90) {
    for (j in 1:90) {
        if (mymat[i,j] >0) {
            ni <- c(ni, rownames(mymat)[i])
            nj <- c(nj, colnames(mymat)[j])

        }

    }
}
list(ni,nj)
}

withab_tmp <- helper_edges(spiec.res.a0.withab.cm)
woutab_tmp <- helper_edges(spiec.res.a0.woutab.cm)

length(withab_tmp[[1]])
#[1] 45
nrow(network_vars_min0_withab$original.edglist[, c("Var1", "Var2")])
#[1] 46
# check
filtered_withab_edglist <- network_vars_min0_withab$original.edglist %>% filter(Var1 %in% taxa_common)
nrow(filtered_withab_edglist[, c("Var1", "Var2")])
#[1] 45 OK.

length(woutab_tmp[[1]])
#[1] 44
nrow(network_vars_min0_woutab$original.edglist[, c("Var1", "Var2")])
#[1] 44

# so how many/ which links are in common

common_a0 <- as.vector(spiec.res.a0.withab.cm) * as.vector(spiec.res.a0.woutab.cm)

as.matrix(getOptCov(spieceasi.res.th001_withab))
as.matrix(getOptCov(spieceasi.res.th001_woutab))

# not that interesting after all

## Mortality associations

# lataa =paste0(datadir, "/full_subnet_analysis_figures_res", analysis_date, "_", analysis_ident, ".rds") ja vertaile

res_tmp_withab <- readRDS("output_2020-04-13_nop_withab/raw_output/full_subnet_analysis_figures_res2020-04-13_noplasmid_withantibiot.rds")

res_tmp_woutab <- readRDS("output_2020-04-13_nop_norxj01/raw_output/full_subnet_analysis_figures_res2020-04-13_noplasmid_norxj01.rds")

res_tmp_woutab %>% filter(componentID==4)
res_tmp_withab %>% filter(componentID==4)

otus_sn4 <- (res_tmp_withab$n.components %>% filter(componentID==4))[,"OTUunderscored"]
otus_sn4b <- (res_tmp_woutab$n.components %>% filter(componentID==4))[,"OTUunderscored"] # same

woutab_sn4 <- res_tmp_woutab$cox_linear %>% filter(predictor %in% otus_sn4)
withab_sn4 <- res_tmp_withab$cox_linear %>% filter(predictor %in% otus_sn4)