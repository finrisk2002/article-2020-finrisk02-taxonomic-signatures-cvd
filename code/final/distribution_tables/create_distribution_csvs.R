library(biomformat)
library(microbiome)
library(phyloseq)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(magrittr)

# helper function
write.extra.csv <- function(df, fn, extra_header, ...) {
  data <- file(fn, open='wt')
  on.exit(close(data))
  writeLines(extra_header, con=data)
  write.table(df, data, ...)
}

tbldir <- outdir

pseq.g.d1<- readRDS(genus_level_phyloseq_path)
pseq.p.d1<- readRDS(phylum_level_phyloseq_path)


pseq.p.d1.rel <- transform(pseq.p.d1, "compositional")

pseq.g.d1.rel <- transform(pseq.g.d1, "compositional")


## table 1: phylum level distribution


get_table_columns <- function(pseq_abs, pseq_rel) {

  pseq_rel_core1 <- core(pseq_rel, prevalence=1/100, detection=0.1/100)

  t1.c1 <- taxa_sums(pseq_abs)

  # col 2: mean relative abundance

  t1.c2 <- as.numeric(format(apply(abundances(pseq_rel), 1, mean),
                            width=10, signif=9, scientific=6))

  t1.c2 <- (100*(abundances(pseq_rel) %>% rowMeans())) %>% round(9) %>% format(width=10, signif=9, scientific=FALSE, nsmall=8)

  # col 3: 95% relative abundance interval
  t1.c3 <- apply(abundances(pseq_rel), 1, function(genus) {
    q <- (100*quantile(genus, probs = c(0.025, 0.975))) %>% round(9) %>% format(width=10, signif=9, scientific=FALSE, nsmall=8)

    paste0(q, collapse = "-")
  })


  # col 4: prevalence at detection 0.0

  t1.c4 <- as.numeric(100*prevalence(pseq_rel)) %>% round(9) %>% format(width=10, signif=9, scientific=FALSE, nsmall=8)

  # col 5: prevalence at detection=0.1/100

  t1.c5 <- as.numeric(100*prevalence(pseq_rel, 0.1/100)) %>% round(9) %>% format(width=10, signif=9, scientific=FALSE, nsmall=8)

  t1.c6 <- sapply(names(t1.c1), function (tn) { ifelse(tn %in% taxa_names(pseq_rel_core1),
                                                     "Yes", "No")})

  list(total_read_count=t1.c1,
       mean_relative_abundance=t1.c2,
       sd_relative_abundance=t1.c3,
       prevalence_0det = t1.c4,
       prevalence_01pc_det = t1.c5,
       presence_in_core = t1.c6)
}

# col 1: total read count per phylum

t1_list <- get_table_columns(pseq.p.d1, pseq.p.d1.rel)


t1 <- data.frame(Phylum=names(t1_list$total_read_count),
                 Mean.Relative.Abundance=t1_list$mean_relative_abundance,
                 SD.Relative.Abundance=t1_list$sd_relative_abundance,
                 Prevalence.at.0.1.percent.Detection=t1_list$prevalence_01pc_det)

cat("number of phyla\n")
cat(nrow(t1))
cat("\n")

# Adjust column names and clean phylum names
t1 <- t1 %>%
  set_colnames(c("Phylum",
                 "Mean Relative Abundance",
                 "95% Relative Abundance Interval",
                 "Prevalence at 0.1 Percent Detection"))


## Build and write

write.extra.csv(t1, paste0(tbldir, 'table_shogun_phylum_distribution.csv'),
                 "\"Distribution of phyla in FINRISK 2002\",,,,,\n", sep=",", row.names=FALSE)

## table 2: genus-level + core

t2_list <- get_table_columns(pseq.g.d1, pseq.g.d1.rel)

t2 <- data.frame(Genus=names(t2_list$total_read_count),
                 Mean.Relative.Abundance=t2_list$mean_relative_abundance,
                 Interval.95.Relative.Abundance=t2_list$sd_relative_abundance,
                 Prevalence.at.0.1.percent.Detection=t2_list$prevalence_01pc_det,
                 Present.in.1.percent.Core=t2_list$presence_in_core)
t2old <- t2
t2 <- t2 %>%
  set_colnames(c("Genus",
                 "Mean Relative Abundance (%)",
                 "95% Relative Abundance Interval",
                 "Prevalence at 0.1 Percent Detection (%)",
                 "Included in Prevalence > 1 % Analyses"))

cat("number of genera\n")
cat(nrow(t2))
cat("\n in core \n")
cat(sum(t2_list$presence_in_core == "Yes"))
cat("\n median relative abudnace of core means \n")
cat(median(t2old %>% filter(Present.in.1.percent.Core == "Yes")
  %>% pull(Mean.Relative.Abundance) %>% as.character %>% as.numeric))
cat("\n")
cat("\n combined core median rel abundance \n")

pseq.g.d1.rel.core <- pseq.g.d1.rel %>% core(prevalence=1/100, detection=0.1/100)

cat(median(sample_sums(pseq.g.d1.rel.core)))
cat("\n")


write.extra.csv(t2, paste0(tbldir, 'table_shogun_genus_distribution.csv'),
                 paste0("\"Distribution of genera in FINRISK 2002\",,,,\n"),
                sep=",", row.names=FALSE)

# no species thus far



