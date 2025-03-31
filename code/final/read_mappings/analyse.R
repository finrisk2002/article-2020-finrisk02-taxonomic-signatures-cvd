system(paste("rm coverage.tab"))
path <- "/home/lemila/data/tmp/BAM/"
coverage.files <- list.files(path, pattern = ".coverage$", full.names = TRUE)
hits <- list()

TYHJAT hitit ei näytä tulleen (hits <- c()) -> testaa ja muokkaa ne mukaan
# fs <- sample(coverage.files, 1000)
fs <- coverage.files
for (f in fs) {

  print(paste("analyse.R", match(f, fs)/length(fs)))

  if (length(readLines(f)) > 0) {
    x <- read.table(f)
    hits[[f]] <- x
  } else {
    hits[[f]] <- c()
  }
  if (match(f, fs) %in% seq(500, 1e4, 500)) {
    gc()
  }
}

names(hits) <- gsub(".bam.coverage", "", gsub(path, "", names(hits)))

saveRDS(hits, file = "hits.Rds")



# hits <- lapply(names(hits), function (i) {x <- hits[[i]]; x$Sample <- rep(i, nrow(x))})
# x$Sample <- rep(, nrow(x))    
# names(x) <- c("Gene", "Location", "Depth")

