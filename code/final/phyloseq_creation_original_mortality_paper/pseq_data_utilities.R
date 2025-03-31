# functions for reading and preprocessing various data files and objects for pseq creation
# Aaro Salosensaari May 2018

#' Read sample ('phenotype') data.
#'
#' Reads the phenotype TSV file with fread and sets sample_names
#' to Barcode.
#'
#' @importFrom data.table fread
#' @importFrom phyloseq sample_data sample_names sample_names<-
#'
#' @export
read_pheno_data <- function (pheno_file, stool.id='Barcode') {

  # get phenotype / sample data from file
  pheno_data_raw <- as.data.frame(fread(pheno_file, header=T))

  # skip the rows where the ID (Barcode in 02, STOOL_ID in 07 data) is missing
  pheno_data_na <- is.na(pheno_data_raw[,stool.id])
  pheno_data <- sample_data(pheno_data_raw[!pheno_data_na, -1])
  # sample_names given by Barcode
  sample_names(pheno_data) <- pheno_data_raw[!pheno_data_na, stool.id]
  pheno_data
}


#' Read updated sample ('phenotype') data.
#'
#' Reads the updated phenotype Rdata file with load and sets sample_names
#' to Barcode.
#' (The Rdata file is /csc/fr_metagenome/pheno/2015_60_Salomaa_Jain_dataFR02_2018-07-04.RData as of 2018-08-14.)
#'
#' @importFrom phyloseq sample_data sample_names sample_names<-
#'
#' @export
read_pheno_data2 <- function (pheno_rdata_file, stool.id='Barcode', dataset.name="FR02") {

  # get phenotype / sample data from file
  tmpenv <- new.env()
  load(pheno_rdata_file, tmpenv)

  pheno_data_raw <- tmpenv[[dataset.name]]

  # if Rdata file contains these columns, ignore them
  ignored.cols <- c("WESID", "WGSID", "BATCH", "FID", "DNAPERM")

  pheno_data_raw <- pheno_data_raw[, !colnames(pheno_data_raw) %in% ignored.cols]

  if ("NA" %in% colnames(pheno_data_raw) && !("NA." %in% colnames(pheno_data_raw))) {
    pheno_data_raw[,'NA.'] <- pheno_data_raw[,'NA']
    pheno_data_raw[,'NA'] <- NULL
  }

  # skip the rows where the ID (Barcode in 02, STOOL_ID in 07 data) is missing
  pheno_data_na <- is.na(pheno_data_raw[,stool.id])
  pheno_data <- sample_data(pheno_data_raw[!pheno_data_na, -1])
  # sample_names given by Barcode
  sample_names(pheno_data) <- pheno_data_raw[!pheno_data_na, stool.id]
  pheno_data
}


#' Clean sample ('phenotype') data.
#'
#' * Cleans CRP, CDT to numeric.
#' * Takes dataframe names.df that contains metadata descriptions of the sample
#' covariates, and transforms the variables to numeric/ordered/factor as per the
#' metadata descriptions. Variables wihtout description in names.df are left as-is.
#' * Also adds labels to couple of factor variables where they make sense.
#'
#' @importFrom microbiome meta
#' @importFrom dplyr %>%
#' @importFrom dplyr n_distinct
#'
#' @export
clean_pheno_data <- function (pheno_data, names.df=NULL) {
  # full doc todo
  # 'data/cv_descriptions.csv'



  pheno_data_tmp <- pheno_data

  # NOTE: for some reason default read methods do nophyloseq::t parse CRP, CDT as numeric,
  # and we have to fix it here
  pheno_data_tmp$CRP <- sapply(pheno_data_tmp$CRP, function(x) as.numeric(x))
  pheno_data_tmp$CDT <- sapply(pheno_data_tmp$CDT, function(x) as.numeric(x))


  pheno_data_tmp <- meta(pheno_data_tmp)

  # variables marked as ordinals, categorical in csv to ordered, factors
  if (!is.null(names.df)) {
    for (vi in 1:nrow(names.df)) {
      v.covar <- names.df$Covariate[vi] %>% as.character
      v.type <- names.df$Type[vi] %>% as.character
      cat(v.covar)
      cat('\n')
      cat(v.type)
      cat('\n')
      if (v.type == "N") {
        pheno_data_tmp[, v.covar] <- sapply(pheno_data_tmp[, v.covar],
                                          function (x) as.numeric(x))
      }
      else if (v.type == "O") {
        pheno_data_tmp[, v.covar] <- ordered(pheno_data_tmp[, v.covar])
      }
      else if (v.type == "No") {
        pheno_data_tmp[, v.covar] <- factor(pheno_data_tmp[, v.covar])
      }
      else if (v.type == "D") {
        pheno_data_tmp[, v.covar] <- factor(pheno_data_tmp[, v.covar])
      }
      else {
        # unknown type, assume numeric
        cat(paste("Casting unknown variable type", v.type, "to numeric\n"))
        pheno_data_tmp[, v.covar] <- sapply(pheno_data_tmp[, v.covar],
                                          function (x) as.numeric(x))
      }
    }
  }

  pheno_data_clean <- sample_data(pheno_data_tmp)



  # edit labels of some factors to more descriptive ones
  pheno_data_clean <- map_num_to_named_factors(pheno_data_clean)

  pheno_data_clean
}

#' Process biom files from SHOGUN to phyloseq objects.
#'
#' Additional parameter readable_rownames control how OTU (row)names are
#' set. If TRUE, uses naming scheme "Name (Kingdom)". If FALSE, "sp<rowNO>" is used.
#' Needed for processing profile-level biom without introducing duplicate rownames.
#'
#' If parameter skip_profile_nan is TRUE, drop rows with name "nan". (Also needed
#' for profile-level biom.)
#'
#' @importFrom phyloseq phyloseq taxa_names taxa_names<-
#'
#' @examples
#'
#' # Assumes that the script is run on Atlas and files are available at specified locations.
#' # Change 'genus' for 'phylum','species', etc to produce a phylum level objects
#' # from the respective .biom file.
#' # for full script  see inst/examples/genusbiom2phyloseq.R
#' \dontrun{
#' taxo_dataraw_dir <- "/csc/fr_metagenome/microbiome2018-02-26/taxonomy/"
#' genus_counts_file <- paste0(taxo_dataraw_dir, "shogun/combined_redist.genus.biom")
#' #
#' pheno_file <- "/csc/fr_metagenome/thl_fr2_2016/FR02.tsv" # i.e. sample_data
#' #
#' outdir <- "./data_private/"
#' dir.create(outdir)
#' # descriptions etc metadata for our covariates
#' desc.file <- system.file("extdata", "cv_descriptions.csv", package = "finriskmetagcommon", mustWork = TRUE)
#' names.df <- read_metadata_descriptions(desc.file)
#'
#' # read the phenotype data
#' pheno_data <- read_pheno_data(pheno_file)
#' # transform categorical variables to factors and other post-processing
#' pheno_data_clean <- clean_pheno_data(pheno_data, names.df)
#' #
#' # create Phyloseq-object from biom files and cleaned data
#' phfinrisk_genus_clean    <- process_shogun_biom(genus_counts_file, pheno_data_clean)
#' #
#' # save the object
#' saveRDS(phfinrisk_genus_clean, file = paste0(outdir,"phfinrisk_genus_clean.RDs"))
#' # for future reference, also save names.df and with a more descriptive variable name
#' phfinrisk_metadatadesc <- names.df
#' save(phfinrisk_metadatadesc, file=paste0(outdir, "phfinrisk_metadatadesc.RData"))
#'
#' # this metadata description df is also available in the package for easier use
#' # attach(phfinrisk_metadatadesc)
#'
#' }
#'
#' @export
process_shogun_biom <- function (counts_file, pheno_data, skip_profile_nan=FALSE, readable_rownames=TRUE) {

  otu <- extract_otu_table_biom(counts_file, pheno_data, skip_profile_nan)


  my_tt <- extract_taxtable(otu, readable_rownames)
  taxa_names(otu) <- taxa_names(my_tt)

  phfinrisk <- phyloseq(otu, my_tt, pheno_data)

}


#' Extract otu table from biom
#'
#' Helper for process_shogun_biom
#'
#' @importFrom phyloseq otu_table
#'
#' @export
extract_otu_table_biom <- function (counts_file, pheno_data, skip_profile_nan=FALSE) {
  # read the otu data to phyloseq OTU table format

  counts_raw <- biomformat::read_biom(counts_file)
  counts_data = biomformat::biom_data(counts_raw)

  # Warning message:
  # In strsplit(msg, "\n") : input string 1 is invalid in this locale
  # idk what to do about it

  otu <- otu_table(as(counts_data, 'matrix'), taxa_are_rows=T)

  # parse taxa names

  # combined_profile.biom is special, as it contains an OTU named 'nan'.
  # we usually want to skip it

  if (skip_profile_nan) {
    nonnan_otus <- rownames(otu) != 'nan'
    otu <- otu[nonnan_otus,]
  }

  otu
}

#' Extract tax table from otu names
#'
#' Helper for process_shogun_biom
#'
#' @importFrom phyloseq build_tax_table taxa_names parse_taxonomy_qiime taxa_names<-
#'
#' @export
extract_taxtable <- function (otu, readable_rownames=TRUE) {

  parsed_taxa_names <- lapply(rownames(otu),
                              function(x) parse_taxonomy_qiime(x));

  my_tt <- build_tax_table(parsed_taxa_names)

  if (readable_rownames) {
    # set taxa names to something that can be easily interpreted instead of sp(number)
    # apparently this is no longer recommended but makes usage of microbiome plotting functions easier
    nonnas <- !is.na(my_tt)
    last_nonnas <- apply(nonnas, 1, function (r) {max(which(r))})

    tt_n <- nrow(my_tt)
    tt_names <- sapply(1:tt_n, function (i) {paste0(as.character(my_tt[i,last_nonnas[i]]), " (", as.character(my_tt[i,'Kingdom']), ")")})

    taxa_names(my_tt) <- tt_names
  }

  my_tt
}


#' Read metadata descriptions from csv ( ;-separated file).
#'
#' @importFrom dplyr %>%
#'
#' @export
read_metadata_descriptions <- function (metadata_file,
                                        separator=';') {

  variable.desc <- read_metadata_raw(metadata_file, separator=separator)

  raw.colnames <- colnames(variable.desc)

  extra.colnames <- raw.colnames[!raw.colnames %in% c('VARIABLE',
                                                      'Variable.Name',
                                                      'Variable.Description',
                                                      'Category',
                                                      'Type')]

  # process compulsory columns and create a df

  var.covar <- sapply(1:length(variable.desc$VARIABLE), function (vi) {
                        v.name <- variable.desc$VARIABLE[vi] %>% as.character })

  var.names <- sapply(1:length(variable.desc$VARIABLE), function (vi) {
                        v.name <- variable.desc$Variable.Name[vi] %>% as.character })

  var.longnames <- sapply(1:length(variable.desc$VARIABLE), function (vi) {
                          v.longname <- variable.desc$Variable.Description[vi] %>% as.character
                        })

  var.category <- sapply(1:length(variable.desc$VARIABLE), function (vi) {
                        v.cat <- variable.desc$Category[vi] %>% as.character })

  var.types <- sapply(1:length(variable.desc$VARIABLE), function (vi) {
                        v.cat <- variable.desc$Type[vi] %>% as.character })

  names.df <- data.frame(Covariate=var.covar, Name=var.names,
                         Desc=var.longnames, Category=var.category,
                         Type=var.types)

  # add arbitrary columns as characters
  for (cn in extra.colnames) {
    names.df[,cn] <- sapply(1:length(variable.desc$VARIABLE), function(vi) {
                          v.cn <- variable.desc[vi,cn] %>% as.character
                         })
  }

  # return

  names.df
}

read_metadata_raw <- function (metadata_file, separator=";") {

  # TODO: use fread here too for consistency
  variable.desc <- utils::read.table(metadata_file, header=T, na.strings = "na",
                              strip.white=T, sep=separator, quote="")
  variable.desc
}



#' Recode named factors back to numeric.
#'
#' @param pheno_data phyloseq sample data object
#' @return pheno_data phyloseq sample data object
#'
#' @importFrom dplyr recode
#' @export
map_named_factors_to_num <- function(pheno_data) {

  pheno_data$MEN    <- recode(pheno_data$MEN,
                              "Female"=0, "Male"=1)
  pheno_data$EAST   <- recode(pheno_data$EAST,
                              "WEST"=0,"EAST"=1)
  if ('ALUE' %in% colnames(pheno_data) && n_distinct(pheno_data$ALUE) > 2) {
    pheno_data$ALUE <- recode(pheno_data$ALUE,
                              "North Karelia"=2,
                              "North Savonia"=3,
                              "Turku/Loimaa"=4,
                              "Helsinki/Vantaa"=5,
                              "Oulu province"= 6,
                              "Lapland"=7)
  }
  if ('KOULGR' %in% colnames(pheno_data)) {
    pheno_data$KOULGR <- recode(pheno_data$KOULGR,
                              "low"=1, "medium"=2, "high"=3)
  }
  pheno_data
}


# see map_named_factors_to_num
map_MEN_EAST_factors_to_num <- function(pheno_data) {

  pheno_data$MEN    <- recode(pheno_data$MEN,
                              "Female"=0, "Male"=1)
  pheno_data$EAST   <- recode(pheno_data$EAST,
                              "WEST"=0,"EAST"=1)
  pheno_data
}


#' Recode certain factors that have easily interpretable names.
#'
#' NOTE. Will produce nonsensical results if sample_data has
#' been e.g. scaled
#'
#' @param pheno_data phyloseq sample data object
#' @return pheno_data phyloseq sample data object
#'
#' @importFrom dplyr recode_factor
#' @export
map_num_to_named_factors <- function(pheno_data) {
  pheno_data$MEN    <- recode_factor(factor(pheno_data$MEN),
                                     "0"="Female", "1"="Male")
  pheno_data$EAST   <- recode_factor(factor(pheno_data$EAST),
                                     "0"="WEST","1"="EAST")
  if ('ALUE' %in% colnames(pheno_data) && n_distinct(pheno_data$ALUE) > 2) {
    pheno_data$ALUE <- recode_factor(factor(pheno_data$ALUE),
                                     "2"="North Karelia",
                                     "3"="North Savonia",
                                     "4"="Turku/Loimaa",
                                     "5"="Helsinki/Vantaa",
                                     "6"="Oulu province",
                                     "7"="Lapland")
  }
  if ('KOULGR' %in% colnames(pheno_data)) {
    pheno_data$KOULGR <- recode_factor(factor(pheno_data$KOULGR),
                                     "1"="low", "2"="medium", "3"="high")
  }
  pheno_data
}

#' Ensures all are numeric and scales them.
#'
#' i.e. substracts the mean and divides by standard deviation.
#'
#' @param pseq phyloseq object
#' @param names.df unused argument
#' @return a phyloseq object where sample_data has been scaled
#'
#' @importFrom microbiome meta
#' @importFrom phyloseq sample_data
#'
#' @export
scale_sample_data_to_num <- function (pseq, names.df=NULL) {

  # center and scale the phenotype measurements (by 1sd)
  # all phenotype data variables have now same units and beta is more easier to interpret

  pheno_data <- map_named_factors_to_num(sample_data(pseq))

  vnames <- colnames(meta(pheno_data))
  X <- data.frame(lapply(meta(pheno_data[,vnames]),
                        function (x) {
                          if (is.factor(x)) {
                            scale(as.numeric(levels(x))[x])
                          } else {
                            scale(as.numeric(x))
                          }
                        }))


  rownames(X) <- rownames(pheno_data)

  sample_data(pseq) <- sample_data(X)

  pseq
}




#' Set all variables that are marked dichotomous (type == "D") to factors.
#'
#' @param pseq phyloseq object with numeric sample data
#' @param names.df metadata data.frame
#'
#' @importFrom dplyr filter
#'
#' @export
set_dichotomous_sample_data_to_fact <- function (pseq, names.df) {

  vnames.dicho <- as.character(filter(names.df, Type=="D")[['Covariate']])

  pheno_data <- map_num_to_named_factors(sample_data(pseq))
  vnames_all <- colnames(meta(pheno_data))
  X <- data.frame(lapply(vnames_all,
                        function (vn) {
                          if (vn %in% vnames.dicho) {
                            x <- meta(pheno_data)[,vn]
                            cat(paste("setting", vn, " to factor\n"))
                            factor(x)
                          } else {
                            meta(sample_data(pseq))[,vn]
                          }
                        }))

  colnames(X) <- vnames_all


  rownames(X) <- rownames(pheno_data)

  sample_data(pseq) <- sample_data(X)

  pseq

}


#' Set *non-dichotomous* variables to numeric and scales them.
#'
#' i.e. substracts the mean and divides by standard deviation;
#' dichotomous factor are not touched
#'
#' @param pseq phyloseq object with named levels
#' @param names.df metadata data.frame
#'
#' @return a phyloseq object where sample_data has been scaled
#'
#' @importFrom microbiome meta
#' @importFrom phyloseq sample_data
#'
#' @export
scale_nondicho_sample_data_to_num <- function (pseq, names.df) {

  # center and scale the phenotype measurements (by 1sd)
  # all phenotype data variables have now same units and beta is more easier to interpret

  vnames.dicho <- as.character(filter(names.df, Type=="D")[['Covariate']])

  vnames_all <- colnames(meta(sample_data(pseq)))

  pheno_data_tmp <- map_named_factors_to_num(sample_data(pseq))
  X <- data.frame(lapply(vnames_all,
                        function (vn) {
                          x <- meta(sample_data(pseq))[,vn]
                          if (!vn %in% vnames.dicho) {
                            cat(paste("scaling", vn, "\n"))
                            if (is.ordered(x)) {
                              scale(as.numeric(x))
                            } else if (is.factor(x)) {
                              xt <- meta(pheno_data_tmp)[,vn]
                              scale(as.numeric(xt))
                            } else {
                              scale(as.numeric(x))
                            }
                          } else {
                            x
                          }
                        }))

  colnames(X) <- vnames_all


  rownames(X) <- rownames(sample_data(pseq))

  sample_data(pseq) <- sample_data(X)

  pseq
}

#' Calls scale_nondicho_sample_data_to_num and set_dichotomous_sample_data_to_fact
#'
#' @export
scale_nondicho_set_dicho <- function (pseq, names.df) {
  pseq <- scale_nondicho_sample_data_to_num(pseq, names.df)
  pseq <- set_dichotomous_sample_data_to_fact(pseq, names.df)
  pseq
}



#' Ensures all data are numeric but does not scale them.
#'
#'
#' @param pseq phyloseq object
#' @param names.df unused argument
#' @return a phyloseq object where sample_data has been scaled
#'
#' @importFrom microbiome meta
#' @importFrom phyloseq sample_data
#'
#' @export
set_sample_data_to_num <- function (pseq, names.df=NULL) {

  # center and scale the phenotype measurements (by 1sd)
  # all phenotype data variables have now same units and beta is more easier to interpret

  pheno_data <- map_named_factors_to_num(sample_data(pseq))

  vnames <- colnames(meta(pheno_data))
  X <- data.frame(lapply(meta(pheno_data[,vnames]),
                        function (x) {
                          if (is.factor(x)) {
                            as.numeric(levels(x))[x]
                          } else {
                            as.numeric(x)
                          }
                        }))


  rownames(X) <- rownames(pheno_data)

  sample_data(pseq) <- sample_data(X)

  pseq
}


#' Set *non-dichotomous* variables to numeric but does not scale them.
#'
#' dichotomous factor are not touched
#'
#' @param pseq phyloseq object with named levels
#' @param names.df metadata data.frame
#'
#' @return a phyloseq object where sample_data has been scaled
#'
#' @importFrom microbiome meta
#' @importFrom phyloseq sample_data
#'
#' @export
set_nondicho_sample_data_to_num <- function (pseq, names.df) {

  # center and scale the phenotype measurements (by 1sd)
  # all phenotype data variables have now same units and beta is more easier to interpret

  vnames.dicho <- as.character(filter(names.df, Type=="D")[['Covariate']])

  vnames_all <- colnames(meta(sample_data(pseq)))

  pheno_data_tmp <- map_named_factors_to_num(sample_data(pseq))
  X <- data.frame(lapply(vnames_all,
                        function (vn) {
                          x <- meta(sample_data(pseq))[,vn]
                          if (!vn %in% vnames.dicho) {
                            cat(paste("setting", vn, " to numeric\n"))
                            if (is.ordered(x)) {
                              as.numeric(x)
                            } else if (is.factor(x)) {
                              xt <- meta(pheno_data_tmp)[,vn]
                              as.numeric(xt)
                            } else {
                              as.numeric(x)
                            }
                          } else {
                            x
                          }
                        }))

  colnames(X) <- vnames_all


  rownames(X) <- rownames(sample_data(pseq))

  sample_data(pseq) <- sample_data(X)

  pseq
}


#' Sets non-dichotomous variables to numeric and dichotomous to factors
#'
#' @export
set_nondicho_set_dicho <- function (pseq, names.df) {
  pseq <- set_nondicho_sample_data_to_num(pseq, names.df)
  pseq <- set_dichotomous_sample_data_to_fact(pseq, names.df)
  pseq
}



#' Drop variables in phyloseq sample data that are not in names.df
#'
#' @param pseq_level_all phyloseq object with all possible sample variables
#' @param names.df metadata data file (such as phfinrisk_metadatadesc)
#'
#' @return a phyloseq object, sample variables without row in names.df removed
#'
#' @importFrom phyloseq sample_data
#' @importFrom microbiome meta
#'
#' @export
drop_unused_samplevars <- function(pseq_level_all, names.df) {
  pseq_level_ok <- pseq_level_all
  tmp_samples <- meta(sample_data(pseq_level_all))
  tmp_samples_ok <- tmp_samples[colnames(tmp_samples) %in% names.df$Covariate]
  sample_data(pseq_level_ok) <- sample_data(tmp_samples_ok)
  pseq_level_ok
}


#' Drop samples with small readcounts from a phyloseq.
#'
#' Convenience wrapper for prune_samples.
#'
#' @param pseq phyloseq object
#' @param countthr read count threshold
#'
#' @return phyloseq object with samples with reads <= countthr dropped
#'
#' @importFrom phyloseq sample_sums
#' @importFrom phyloseq prune_samples
#'
#' @export
drop_small_rc_pseq <- function(pseq, countthr=50000) {
  included.samples <- which(sample_sums(pseq) > countthr)
  pseq.dropped <- prune_samples(names(included.samples), pseq)
  pseq.dropped
}

#' Rarefy and prune a phyloseq at specified depth.
#'
#' Wrapper for phyloseq::rarefy_even_depth that also drops
#' the samples fall below the specified read count threshold.
#'
#' @param pseq phyloseq object
#' @param countthr read count threshold
#' @param rngseed rngseed, passed to rarefy_even_depth
#'
#' @return phyloseq object with samples with reads <= countthr dropped  and
#'         rest rarefied
#'
#' @importFrom phyloseq rarefy_even_depth
#'
#' @export
rarefy_pseq <- function(pseq, countthr=50000, rngseed=42) {

  pseq.dropped <- drop_small_rc_pseq(pseq, countthr=countthr)
  pseq.rarefied <- rarefy_even_depth(pseq.dropped,
                                     sample.size=countthr,
                                     rngseed=rngseed)


  pseq.rarefied
}


#' Creates a Phyloseq otu_table from Centrifuge bioms
#'
#' The said bioms being  /csc/fr_metagenome/microbiome2018-02-26/taxonomy/centrifuge/combined*biom
#' Additionally can attempt to parse some information from the kraken reports directly (slow)
#'
#' @param counts_file biom file that contains abundance counts
#' @param reportdirs_named named vector of dirs that contain reports from which pretty 
#'        human-readable taxa names can be extracted. names(reportdirs_named) should be
#'        sample ids.
#' @param extract_nice_names logical, whether to use reportdirs_named to get the nice
#'        human-readable names
#'
#' @export
#'
#' @importFrom hashmap hashmap
#' @importFrom phyloseq otu_table
extract_otu_table_centrifuge_biom <- function(counts_file, reportdirs_named=NULL, extract_nice_names=FALSE) {

  counts_raw <- biomformat::read_biom(counts_file)
  counts_data_old = biomformat::biom_data(counts_raw)
  counts_data <- counts_data_old

  if (!is.null(reportdirs_named) || extract_nice_names) {
    snames <- colnames(counts_data)

    H <- hashmap('_','')
    H$clear()

    for (si in 1:length(snames)) {
      sn <- snames[si]
      dpath <- reportdirs_named[[sn]]
      cat(dpath)
      cat(" ... ")
      cat(si)
      cat("/")
      cat(length(snames))
      cat(" ...\n")
      txt <- read.table(file.path(dpath, "centrifuge", paste0(sn, ".report.txt")), sep='\t')
      # TODO: we could also check other content of report / profile.txts
      taxIDs <- txt[, 5]
      unclean_names <- txt[, 6]
      clean_names <- sapply(unclean_names, trimws)
      already_mapped <- H$has_keys(taxIDs)
      H[[taxIDs[!already_mapped]]] <- clean_names[!already_mapped]
    }
    rownames(counts_data) <- H[[rownames(counts_data_old)]]
  }

  otu <- otu_table(as(counts_data, 'matrix'), taxa_are_rows=T)
}


#' Process biom files from Centrifuge to phyloseq objects.
#'
#' The said bioms being  /csc/fr_metagenome/microbiome2018-02-26/taxonomy/centrifuge/combined*biom
#' TODO: parse all information from the kraken reports directly
#'
#' @param counts_file biom file that contains abundance counts
#' @param pheno_data phenotype sample data (will be set as sample_data(pseq))
#' @param reportdirs_named named vector of dirs that contain reports from which pretty 
#'        human-readable taxa names can be extracted. names(reportdirs_named) should be
#'        sample ids.
#' @param extract_nice_names logical, whether to use reportdirs_named to get the nice
#'        human-readable names
#'
#' @seealso extract_otu_table_centrifuge_biom process_shogun_biom
#'
#' @importFrom phyloseq phyloseq taxa_names taxa_names<-
#'
#' @export
process_centrifuge_biom <- function (counts_file, pheno_data, reportdirs_named=NULL, extract_nice_names=FALSE) {

  otu <- extract_otu_table_centrifuge_biom(counts_file, reportdirs_named=reportdirs_named, extract_nice_names=extract_nice_names)
  #  NOTE: tax table could be parsed from reports, but we will now skip this
  #my_tt <- extract_taxtable(otu, readable_rownames)
  #phfinrisk <- phyloseq(otu, my_tt, pheno_data)
  phfinrisk <- phyloseq(otu, pheno_data)
}
