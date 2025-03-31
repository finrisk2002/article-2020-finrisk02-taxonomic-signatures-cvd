# Readme

This subdirectory contains helpful scripts and other information  to create phyloseq objects from the biom-files / Kraken-style SHOGUN/Centrifuge report files provided by Knight Lab.

## Collecting a biom file from Centrigue (kraken-like) report files

For some sets (SHOGUN, Centrifuge), biom was already provided, so the steps underneath are not strictly necessary but recorded if future user finds them useful. (Also note that SHOGUN output is not given in kraken-like report format.)

For this there is several approaches (including some custom R code). Currently the easiest to use is probably kraken-biom: https://github.com/smdabdoub/kraken-biom

### Instructions for installing kraken-biom

There are many ways to do this. One way to do it is as follows:

1. Install miniconda3 https://docs.conda.io/en/latest/miniconda.html
2. In miniconda3, create an environment (see `conda create`, `conda activate`), install biom-format, then install kraken-biom https://github.com/smdabdoub/kraken-biom .

### Running kraken-biom on Centrifuge (phylum to species level)

`kraken-biom --max P --min S /csc/fr_metagenome/microbiome2018-02-26/taxonomy/*/centrifuge/*.report.txt -o centrifuge_table.biom`

## Creating phyloseqs from bioms and/or other Shogun outputs

Master example script for creating phyloseqs from various data sources is `phyloseq_creation.R`,
which calls functions defined in `pseq_data_utilities.R`.

This master script is not strictly necessary either (see https://joey711.github.io/phyloseq/import-data.html), but can be useful for e.g. setting dichotomous data variables to R factors, and setting prettier genus names.

## List of covariates needed for analysis

See also covariate_descriptions.csv

Type for variables listed below can be either Numeric or Dichotomous.

```
VARIABLE    Variable.Name   Variable.Description    Category    Type    LONGNAME
BL_USE_RX_L Cancer&Immuno Rx    Antineoplastic and immunomodulating agents (ATC: L) Medication  D   
CURR_SMOKE  Smoking Current smoker  Lifestyle   D   Current smoker 1=yes
EAST    Region  Living in Eastern vs. Western Finland (0=West, 1=East)  Demographic D   Region, 1=EASTern Finland, 0=WESTern Finland
BP_TREAT    BP Rx   Blood pressure medication at baseline   Medication  D   BP medication at baseline
MEN Sex Sex (0=Women, 1=men)    Demographic D   Men=1, Women=0
DEATH   DEATH   Dead=1  Unknown N   Dead=1
DEATH_AGEDIFF   DEATH_AGEDIFF   Death-age – baseline age  Unknown N   Death-age – baseline age
PREVAL_DIAB Diabetes    History of type I or type II diabetes   Disease D   na
BL_AGE  Age Age at baseline Demographic N   Baseline age
BMI BMI Body mass index Physical    N   Body mass index
SYSTM   Systolic BP Systolic blood pressure Physical    N   Systolic bp
```

also required for some sensitivity analyses:
```
BL_USE_RX_J01   RX_J01  Antibacterials for systemic use  (ATC: J01) Medication  N
HFC_score   HFC_score   Healthy food choices score  Diet    N
```