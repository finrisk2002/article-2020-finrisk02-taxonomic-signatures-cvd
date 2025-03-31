# Source code

# IN PROGRESS

Material for taxonomic analysis of FINRISK 2002 with CVD endpoint.

Based on material from https://gitlab.com/finrisk2002/article_2019_finrisk02_taxonomic_signatures_mortality_risk

Taxonomic Signatures of Long-Term Mortality Risk in Human Gut
Microbiome.  Aaro Salosensaari, Ville Laitinen, Aki S Havulinna,
Guillaume Meric, Susan Cheng, Markus Perola, Liisa Valsta, Georg
Alfthan, Michael Inouye, Jeramie D Watrous, Tao Long, Rodolfo Salido,
Karenina Sanders, Caitriona Brennan, Gregory C Humphrey, Jon G
Sanders, Mohit Jain, Pekka Jousilahti, Veikko Salomaa, Rob Knight, Leo
Lahti, Teemu Niiranen. medRxiv 2019.12.30.19015842; doi:
https://doi.org/10.1101/2019.12.30.19015842

**Contact:** Leo Lahti leo.lahti@utu.fi for any queries regarding the
  analysis code.

## Comment

This repository contains R scripts and functions necessary (+ some extra) for reproducing analyses with the data given in formats defined in setup scripts (see below).

## People

TBA

## System requirements

Prerequisites for statistical analyses: R, recommended version 3.6.3 or later.

### Installation of R package dependencies

Install R https://cran.r-project.org/, and run the following R code which installs any missing packages and outputs their versions:

```
inst_packages <- c("phyloseq", "microbiome", "survival",
                    "biomformat",  "vegan",
                    "dplyr", "tidyverse", "magrittr", "reshape2",
                    "ggplot2", "scales",
                    "network", "sna", "GGally",
                    "RColorBrewer",
                    "tibble",
                    "knitr", "smoothHR", "openxlsx", 
                    "randomForestSRC", "survivalROC", 
                    "ggfortify", "ggsignif",
                    "grid", "gridExtra", "cowplot", "gtable", "ggpubr", "ggsci",
                    "devtools")

missing_packages <- inst_packages[!(inst_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) != 0) install.packages(missing_packages)

# for network analyses, we use SpiecEasi, which can not be installed with install.packages

# see https://github.com/zdk123/SpiecEasi

library(devtools)
if (! ("SpiecEasi" %in% installed.packages()[,"Package"]))
    install_github("zdk123/SpiecEasi")

all_packages <- c(inst_packages, "SpiecEasi")
for (pkg in all_packages){
    print(pkg)
    print(packageVersion(pkg))
}

```
Version number output on our devoloping environment (FIMM Atlas cluster) is recorded in file pkgVersions.txt

NB. Installation time with regular desktop computer may take some time.

## How to run the analyses

The code for the analyses and figures is in the folder [code/final/](code/final/).

All main analyses can be executed with the R code `code/final/main.R`.

Main analyses require access to the FINRISK 2002 data, in particular, the phyloseq objects stored in RDS format. These are distributed separately via THL Biobank.
Download the files to directory `input/data_main/`.

For details, see below.

### Data directories

Scripts expect to find / write things in `./input/data_<...>` and
`./output/raw_output`. The read and write paths are set up accordingly
(by `code/final/main_setup.R`) if sourcing the `code/final/main.F` script.

Main scripts are identified in `code/final/main.R`; ie. sourcing `main.R` will then source the scripts that will run all main survival or network analyses (respectively).

The script `main_setup.R` sources a script (corresponding to various configuration options) which sets up variables pointing to file paths to phyloseqs and other files.
See `main_setup_<...>.R` for detailed list of files required.

### Phyloseq objects

The statistical analyses code requires metagenome counts data and phenotype data (also called metadata) provided as Phyloseq objects ( https://joey711.github.io/phyloseq/ ) in R data format (RDS).

Creation of the Phyloseq objects from the original Shogun output is detailed in directory `code/final/phyloseq_creation/`.

### Reproductibility

The quantitative results of manuscript (figures and tables) can be reproduced running the scripts given the Phyloseq objects. The scripts write them into directory `output`.

Expected run time on regular desktop: Unfortunately PCA, PCoA and Random Forests are quite slow on the real >7k dataset samples, running all analyses can take several hours.

### Logic

The project structure follows the template https://gitlab.com/openresearchlabs/project_template .


### Instruction specific to running code on review.fimm.fi, as of 2020-02-14

The project directory, located on `/data` contains both code and data.

Please copy the materials from `/data` to your home directory (where you should have write access), e.g. like this:

`cd ~`

`cp -r /data/fr02_mortality_review/article_2019_finrisk02_taxonomic_signatures_mortality_risk ./`

Then navigate to the directory

`cd ~/article_2019_finrisk02_taxonomic_signatures_mortality_risk`

and you are all set to follow the instructions below. <(^_._^)>




