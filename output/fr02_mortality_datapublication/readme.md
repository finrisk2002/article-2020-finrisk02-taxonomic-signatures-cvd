# README

## Directory structure:

* `functional/`:
    * contains module, pathway and ko level data (csv)
    * for ko, extra information on metabolite terms (csv and txt; txt is original, csv easier to read in R)
* `phenotype/`:
    * phenotype data (csv)
* `shogun_tables/`:
    * contains (redistributed) counts for 7211 samples of mortality paper for phylum, genus, species levels (csv)
* `phyloseq/`:
    * contains phyloseqs corresponding to genus and species in `shogun_tables/` (in RDS format)
* `qc/`:
    * the original sequencing and qc data products, organized as follows:
        * `qc/fr02` 7231 samples from FINRISK 2002, used in the mortality paper
        * `qc/antibiot_trial` 114 samples from the antibiotics trial 
        * `qc/fr07` 304 samples from FINRISK 2007
        * `qc/technical` technical samples and other data
        * each directory contains directories (named by numeric sample id) that contains subfolders:
            * `(sample id)/atropos_trimmed`
            * `(sample id)/fastqc_per_sample`
            * `(sample id)/filtered`
        * (sample id) refers to 'barcode' in phenotype data.


## Reading csvs with R

If using R, use the following:

`object <- read.csv(file, header=TRUE, stringsAsFactors=TRUE, check.names=FALSE, row.names=1)`

If you have downloaded the source code release, see `code/final/phyloseq_creation/phyloseq_creation_from_csv.R` for creating phyloseqs from CSVs.