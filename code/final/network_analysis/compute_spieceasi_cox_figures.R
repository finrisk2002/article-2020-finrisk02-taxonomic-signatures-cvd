#compute_spieceasi_cox_figures.R 
# runs main_subnet_analysis_figures.R and full_subnet_analysis_figures.R


source("code/final/network_analysis/main_subnet_analysis_figures.R")
source("code/final/network_analysis/full_subnet_analysis_figures.R")
source("code/final/network_analysis/subnet_analysis_figures_onlynet.R")
source("code/final/network_analysis/east_west_subnet_analysis_figures.R")


phfg.pseq.d50k <- readRDS(genus_level_phyloseq_path)

spieceasi.res.th001 <- readRDS(paste0(datadir,
	"/spieceasi_pseq_2018-12-21_drop50k_core_genus_glasso_nl30_th001.RDs"))

# skip the scripts that only are required to produce figures

if (!exists("create_figures") || create_figures == TRUE) {
	tmp <- main_subnet_analysis_figures(phfg.pseq.d50k, datadir, pltdir, spieceasi.res.th001, fdr_levels=1)
	tmp2 <- full_subnet_analysis_figures(phfg.pseq.d50k, datadir, pltdir, spieceasi.res.th001, tmp$pseq_cox_res, fdr_levels=1)
} else {
	tmp2 <- full_subnet_analysis_figures(phfg.pseq.d50k, datadir, pltdir, spieceasi.res.th001, savefigs=FALSE, fdr_levels=1)
}

# obtain n.components and save
saveRDS(tmp2$n.components, file = network_identities_path)

# save tmp2 in general for further use

saveRDS(tmp2, file=paste0(datadir, "/full_subnet_analysis_figures_res", analysis_date, "_", analysis_ident, ".rds"))


if (!exists("create_figures") || create_figures == TRUE) {
	tmp3 <- nocolor_subnet_analysis_figure(phfg.pseq.d50k, datadir, pltdir, spieceasi.res.th001, tmp$pseq_cox_res)
}
# east / west comparison / verification


# 1. obtain east/west specific spieceasi results


spieceasi.res.th001.east <- readRDS(paste0(datadir,
	"/spieceasi_pseq_2018-12-21_drop50k_core_genus_glasso_nl30_th001_east.RDs"))
spieceasi.res.th001.west <- readRDS(paste0(datadir,
	"/spieceasi_pseq_2018-12-21_drop50k_core_genus_glasso_nl30_th001_west.RDs"))



# 2. call east_west_subnet_analysis_figure. it will ...

if (!exists("create_figures") || create_figures == TRUE) {
	tmp_ew <- east_west_subnet_analysis_figure(phfg.pseq.d50k, datadir, pltdir,
		spieceasi.res.th001.east, spieceasi.res.th001.west, fdr_levels=1)
}

# 3. ...split the pseq of "microbial core" to east/west when applicable (Sorry the code is messy)
# 4. ...run the Cox analysis + create basic ggplot objects separately (with full_subnet .... )
# 5. ...then will combine the ggplot objects and save the plot