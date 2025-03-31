# recommended on Atlas to symlink data directories etc to "/csc/fr_metagenome/microbiome_scratch/scratch/data_<you>"

# other directories from the template (not all are used)
analysis_date <- paste0(Sys.Date())
analysis_ident <- "noplasmid_withantibiot"
datadir <- paste0("output/raw_output/", analysis_date, "/", analysis_ident, "/")

# NB. needs to point to directory where main input data files are located
#input_datadir <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/"
input_datadir <- "input/data_main/"



# Paths for the phyloseq objects, saved in .RDS format
phylum_level_phyloseq_path <- paste0(input_datadir,
    "phfinrisk_phylum_ok_drop50k_2018-12-21_nop_hfc.rds")
genus_level_phyloseq_path <- paste0(input_datadir,
    "phfinrisk_genus_all_drop50k_2018-12-21_nop_hfc.rds"
species_level_phyloseq_path <- paste0(input_datadir,
    "phfinrisk_species_all_drop50k_2018-12-21_nop_hfc.rds")
network_identities_path <- paste0(datadir,"/net_components.rds") # if missing, run compute_spieceasi_cox_figures.R
functional_activities_path <- paste0(input_datadir, "function_activity_norxj01.rds")
ko_metabolite_terms_path <- paste0(input_datadir, "/ko_metabolite_terms.txt")

# NB. in some scripts there might be some hard-coded paths that have escaped our attention.
# they should still point to these paths, however.


#input_datadir <- "input/data_work/"

if (!dir.exists(datadir)) {
  dir.create(datadir, recursive=TRUE)
}

#input_dataraw <- "input/data_raw/" # not currently used
#input_datapr <- "input/data_processed/" # not currently used


create_figures <- TRUE # set FALSE because review env does not support X11.

run_functional_cox <- FALSE


