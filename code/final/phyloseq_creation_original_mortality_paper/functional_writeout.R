# write function data out as csvs
# see also phyloseq_writeout.R and survival_analysis/functional_cox.R

library(phyloseq)
library(microbiome)
library(dplyr)
library(readr)
library(magrittr)
library(tidyverse)

functional_A <- readRDS(functional_activities_path)

metabolite_names <- list(modules = data.frame(Metabolite = functional_A[[1]] %>%
                                                rownames,
                                              Alias = paste0("m", 1:nrow(functional_A[[1]]))),
                         pathways = data.frame(Metabolite = functional_A[[2]] %>%
                                                 rownames,
                                               Alias = paste0("p", 1:nrow(functional_A[[2]]))),
                         ko = data.frame(Metabolite = functional_A[[3]] %>%
                                           rownames,
                                         Alias = paste0("k", 1:nrow(functional_A[[3]])))
)


# Set alias colnames
functional_A <- lapply(names(functional_A), function(x) {
  
  df <- functional_A[[x]]
  
  df <- df %>%
    t %>% 
    set_colnames(metabolite_names[[x]]$Alias)
  
  
  return(df)
  
  
}) %>% set_names(names(functional_A))

# get ko group terms
ko_group_terms <- read_file(file = ko_metabolite_terms_path)
ko_group_terms <- ko_group_terms %>%
  str_split("\n") %>% 
  sapply(function(x) str_split(x, "   ")) %>%
  do.call(rbind, .) %>% 
  as_tibble() %>% 
  separate(V1, sep = ":", into = c("KO", "KO_term")) %>% 
  separate(V2, sep = ";", into = c("code", "KO_group")) %>% 
  select(-c(KO, code)) %>% 
  mutate(KO_group = trimws(KO_group))


# extract sample names
pseq_genus <- readRDS(genus_level_phyloseq_path)


# writeout
functional_outdir <- "output/fr02_mortality_datapublication/functional"

if (!dir.exists(functional_outdir)) {
  dir.create(functional_outdir, recursive=TRUE)
}

modules <- data.frame(functional_A$modules)
rownames(modules) <- rownames(meta(pseq_genus))
write.csv(modules, file=paste0(functional_outdir,"/modules.csv"), row.names=TRUE)

pathways <- data.frame(functional_A$pathways)
rownames(pathways) <- rownames(meta(pseq_genus))
write.csv(pathways, file=paste0(functional_outdir,"/pathways.csv"), row.names=TRUE)  

ko <- data.frame(functional_A$ko)
rownames(ko) <- rownames(meta(pseq_genus))
write.csv(ko, file=paste0(functional_outdir,"/ko.csv"), row.names=TRUE)

write.csv(ko_group_terms, file=paste0(functional_outdir,"/ko_metabolite_terms.csv"))