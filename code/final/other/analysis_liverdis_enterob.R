# analyze cross-sectional association between PREVAL_LIVERDIS and Enterobactericea

library(phyloseq)
library(microbiome)
library(dplyr)

if (!exists("analysis_ident"))
    source("code/final/main_setup.R")

source("code/final/R_functions.R")
source("code/final/R_functions/taxa_names_functions.R")

pseq <- readRDS(genus_level_phyloseq_path)

pseq.rel <- transform(pseq, "compositional") 
taxa_names(pseq.rel) <- taxa_names_underscored(pseq.rel)


pseq.rel.clr <- transform(pseq, "compositional") %>% transform("clr")
taxa_names(pseq.rel.clr) <- taxa_names_underscored(pseq.rel.clr)

tmp <- as.data.frame(tax_table(pseq))
rownames(tmp) <- characters_underscored(rownames(tax_table(pseq)))
entero_tt <- tmp[which(tmp$Family == "Enterobacteriaceae"),]
entero_genera <- rownames(entero_tt)


data <- t(abundances(pseq.rel.clr)[entero_genera,])

covariates <- c("BL_AGE","BMI", "MEN",
                "CURR_SMOKE", "PREVAL_DIAB", "BL_USE_RX_L", 
                "SYSTM", "BP_TREAT")

data <- cbind(data, meta(pseq.rel.clr)[, c(covariates, "PREVAL_LIVERDIS")])

cc <- complete.cases(data)
datacc <- data[cc,]
sum(datacc[, "PREVAL_LIVERDIS"])

form <- paste0("PREVAL_LIVERDIS ~ ", paste(c(covariates, entero_genera), collapse = " + "))

# enterobacteria as separate, CLR

# logistic regression. LIVERDIS ~ covariates + enterob.

datacc_s <- datacc
datacc_s[, c("BL_AGE", "BMI", "SYSTM", entero_genera)] <- scale(datacc_s[, c("BL_AGE", "BMI", "SYSTM", entero_genera)])
res <- glm(as.formula(form), data=datacc_s, family=binomial(link="logit"))

# enterobacteria vs all else, two class CLR = constant * logit
# entero ~ covariates + LIVERDIS

abud_entero_sum <- rowSums(t(abundances(pseq.rel)[entero_genera,]))
abud_entero_CLR <- log(abud_entero_sum) - rowMeans(cbind(log(abud_entero_sum), log(1 - abud_entero_sum)))

data2 <- cbind(EnterobacteriaceaeCLR = abud_entero_logit, meta(pseq.rel.clr)[, c(covariates, "PREVAL_LIVERDIS")])

data2_s <- data2
data2_s[, c("BL_AGE", "BMI", "SYSTM", "EnterobacteriaceaeCLR")] <- scale(data2_s[, c("BL_AGE", "BMI", "SYSTM", "EnterobacteriaceaeCLR")])


cc <- complete.cases(data2)
data2cc <- data2[cc,]

data2cc_s <- data2_s[cc,]

form2 <- paste0("EnterobacteriaceaeCLR  ~ ", paste(c(covariates, "PREVAL_LIVERDIS"), collapse = " + "))

res2 <- lm(as.formula(form2), data=data2cc)
res2s <- lm(as.formula(form2), data=data2cc_s)