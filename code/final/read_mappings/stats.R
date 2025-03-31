library(dplyr)
hits <- readRDS("hits.Rds")
names(hits) <- gsub("^/", "", names(hits))

# Coverages per gene
cnt <- 0
stat <- NULL
for (ss in names(hits)) {

  print(paste("stats.R", ss, match(ss, names(hits))/length(hits)))
  
  h <- hits[[ss]]
  names(h) <- c("Gene", "Location", "Depth")  

  genes <- as.character(unique(h$Gene))
  
  for (gene in genes) {

    # print(paste("....", gene))

    # Reads on this gene
    a0 <- subset(h, Gene == gene) %>% arrange(Location)

    # Then remove bases with zero depth to focus on continuous regions
    a <- dplyr::filter(a0, Depth > 0)

    # if (any(diff(a$Location) > 1)) {stop("Split")}

    inds <- c(1, which(diff(a$Location) > 1), nrow(a))

    # Target gene sequence length
    N <- diff(range(a0$Location))         

    # Continuous coverage lengths.
    # Add +1 to the last due to indexing
    covs <- diff(inds)
    covs[[length(covs)]] <- covs[[length(covs)]] + 1
    
    # Number of covered basepairs
    coverage.abs <- covs
    ave.contig.length <- mean(coverage.abs)
    max.contig.length <- max(coverage.abs)    
    N.contig <- length(coverage.abs)

    # Covered fraction of the target sequence
    coverage.rel <- sum(coverage.abs)/N  
    if (!sum(coverage.abs) == nrow(dplyr::filter(a0, Depth > 0))) {stop("Investigate")}

    # Average depth in the covered area
    ave.depth <- mean(subset(a, Depth > 0)$Depth)

    stat <- rbind(stat, c(sample = ss,
                     gene = gene,
		     target.length = N,		     
		     #coverage.abs = coverage.abs,
		     ave.contig.length = ave.contig.length,
		     max.contig.length = max.contig.length,		     
		     N.contig = N.contig,
		     coverage.rel = coverage.rel,
		     ave.depth = ave.depth))

  }

  if (match(ss, names(hits)) %in% seq(2000, 1e4, 2000)) {
    gc()
    saveRDS(stat, file = "tmp.Rds")
  }
  
}

stat2 <- as.data.frame(stat)
for (nam in colnames(stat2)) {
  stat2[[nam]] <- as.character(stat2[[nam]])
}
for (nam in c("target.length","ave.contig.length","max.contig.length","N.contig","coverage.rel","ave.depth")) {
  stat2[[nam]] <- as.numeric(stat2[[nam]])
}

saveRDS(stat2, file = "stat.Rds")


