library(microbiome)

# Load things

dflist <- readRDS("dflist.Rds")
hits <- readRDS("hits.Rds"); names(hits) <- gsub("^/", "", names(hits))
stat <- readRDS("stat.Rds"); stat[, "sample"] <- gsub("^/", "", stat[, "sample"])

# Harmonize names
names(hits) <- gsub("^/", "", names(hits))
stat[, "sample"] <- gsub("^/", "", stat[, "sample"])

# Keep only those samples that are in dflist (original set)
hits <- hits[names(hits) %in% names(dflist)]
stat <- stat[stat[, "sample"] %in% names(dflist), ]

# Phyloseq object for FR02 data
gen <- readRDS("phfinrisk_centrifuge_genus_d50k_2019-06-02.rds")

Norig <- length(dflist)
Nhits <- length(hits)

set.seed(54)
coms0 <- intersect(intersect(names(dflist), names(hits)), stat[, "sample"])
coms <- coms0
# coms <- sample(coms0, 2000)
# coms <- setdiff(coms0, coms2)

dflist <- dflist[coms]
hits <- hits[coms]
stat <- stat[stat[, "sample"] %in% coms,]
gc()

# Number of unique VFDB hits per FR02 sample
l  <- sapply(dflist, function (x) {length(unique(x$rname_header))});
h  <- h2 <- unlist(sapply(dflist, function (df) {df$rname_header}))

print("SSS")

sss <- c()
for (i in 1:length(hits)) {
  print(i)
  sss[[i]] <- length(unique(hits[[i]][,1]))
}
sss <- data.frame(n = sss)
#gc()

print("SSS2")

stat <- as.data.frame(stat)
for (nam in names(stat)[3:7]) {
  stat[, nam] <- as.numeric(as.character(stat[, nam]))
}

print("Enrichment analysis")
source("enrich.R")

print("report")
library(rmarkdown)
library(knitr)
library(ggplot2)
knit("summary.Rmd")
# render("summary.Rmd")


# Calculate prevalence
spl <- split(accepted.hits[, "sample"], accepted.hits[, "gene"])
prev.abs <- sapply(spl, function (x) {length(unique(x))})
prev <- 100 * prev.abs/Norig

# Calculate association with the subnet
subnet <- paste("g_", c("Klebsiella", "Enterobacter", "Lambdavirus", "Salmonella", "Shigella", "Citrobacter", "Escherichia"), sep = "")

# Add tax table
dd <- data.frame(genus = taxa(gen))
tt <- tax_table(dd)
rownames(tt) <- taxa(gen)
gen <- merge_phyloseq(gen, tt)

# Merge subnet and get CLR abundances
#phy <- merge_taxa2(gen, taxa = paste("g_", subnet, sep = ""), name = "subnet")
#a <- abundances(transform(phy, "clr"))["subnet",]
a <- colSums(abundances(microbiome::transform(gen, "compositional"))[subnet,])
sa <- sort(a)

mean.subnet.abundance <- c()
for (gene in unique(accepted.hits[, "gene"])) {
  # Samples for this gene
  samples <- accepted.hits[which(accepted.hits[, "gene"] == gene), "sample"]
  mean.subnet.abundance[[gene]] <- mean(sa[samples])
}


df <- data.frame(Gene = accepted.hits$rname_header[match(names(prev), accepted.hits$gene)],
                 Prevalence = round(prev, 2)
		 ) %>%
      arrange(desc(Prevalence)) 

dff <- df %>% filter(Prevalence > (100 * 10/Norig))
write.table(dff[, c("Prevalence", "Gene")], file = "VFDB_Table.csv", sep = "\t", quote = FALSE, row.names = FALSE)


# covariates <- c("BL_AGE","BMI", "MEN", "CURR_SMOKE", "PREVAL_DIAB", "SYSTM", "BL_USE_RX_L", "BP_TREAT")
