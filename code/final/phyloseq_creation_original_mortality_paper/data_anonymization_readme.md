# FINRISK 2002 gut microbiome demonstration data 

**Contact:**

 * THL/Biobank: admin.biobank (at) thl.fi
 * Maintainer: Leo Lahti leo.lahti (at) utu.fi



## Summary of the data set 

The data set contains summarized taxonomic abundance tables for 6834
FINRISK 2002 samples at the phylum (51), genus (1754), and species
(5304) level. The original sequencing data is not included in this
anonymized data package.

The taxonomic profiling is based on SHOGUN, as described in "Taxonomic
Signatures of Long-Term Mortality Risk in Human Gut Microbiota Aaro
Salosensaari, Ville Laitinen, Aki S Havulinna, Guillaume Meric, Susan
Cheng, Markus Perola, Liisa Valsta, Georg Alfthan, Michael Inouye,
Jeramie D Watrous, Tao Long, Rodolfo Salido, Karenina Sanders,
Caitriona Brennan, Gregory C Humphrey, Jon G Sanders, Mohit Jain,
Pekka Jousilahti, Veikko Salomaa, Rob Knight, Leo Lahti, Teemu
Niiranen" medRxiv 2019.12.30.19015842; doi:

https://doi.org/10.1101/2019.12.30.19015842


### Phenotype data

The data set contains the following phenotype variables for 6834 FINRISK 2002 samples.

 * "EAST": "EAST"/"WEST". Region of subjects' residence (Western or Eastern Finland). 
 
 * "MEN" Sex of the subject. 

 * "BL_AGE_cat" Baseline age categorized into age groups 24_39 (24-39 years), 40_59 (40-59 years) and 60_75 (60-75 years).

 * "BMI_cat": Baseline BMI categorized as "lean" (< 25), "overweight" (>=25 and <30) and "obese" ( >= 30)

 * "DEATH" Whether subject died during followup 15 year follow-up after 2002. 0/1.
 
 * "DEATH_AGEDIFF_anon" Anonymized time (in years) between the baseline time
   point and the last follow-up time point. Random normally
   distributed noise (mean=0, sd=3) was added to the original variable
   "DEATH_AGEDIFF" for the sample that had a reported event, and the
   numbers were rounded to the nearest integer, ensuring that events
   are in range 0 .. 15 . I.E. "DEATH_AGEDIFF_anon" has resolution of
   one year. Effectively, all censoring times without event (original
   range 14.7 .. 14.9) are combined into single integer category (15),
   and each anonymized event time has approximately 31% chance of
   being more than +/- 3 years off from the true quantity.

## Anonymization procedure

Taxonomic profiling data (7211 samples) used in Salosensaari et al. 2020 was
anonymized as follows.

 1. All samples with one or more missing variables were dropped (117 samples).
 
 2. Sample ids were randomly shuffled and renamed s1, s2, etc.
 
 3. Continuous variables (BL_AGE, BMI, DEATH_AGEDIFF) were categorized as described above
 
 4. Discrete variables "EAST", "MEN" and "DEATH" were retained as is.
 
 5. All unique combinations of "EAST", "MEN", "BL_AGE","BMI" and
    "DEATH" with less than 30 samples in each subgroup were removed
    (260 samples). DEATH_AGEDIFF_anon is continuous, but it can be
    characterized as follows: If one looks at the original
    DEATH_AGEDIFF stratified to 5 year intervals, all subgroups have
    at least 6 individuals that are indistinguishable from each other.

## Usage


Both phyloseqs and csv format provided. CSVs can be read in R with command
```
    pheno <- read.csv("phenotype/phenotype.csv")
```

Phyloseqs can be read with

```
    pseq <- readRDS("phyloseq/pseq_genus.rds")
```
