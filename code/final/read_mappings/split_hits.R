library(dplyr)
hits <- readRDS("hits.Rds")
names(hits) <- gsub("^/", "", names(hits))
gc()
w <- 500
inds <- seq(w, length(hits), length = 10)
for (ind in inds) {
  ff <- paste0("subhit", ind, ".Rds", collapse = "")
  print(paste("Writing file", ff))
  subhit <- hits[(ind - w + 1):ind]
  saveRDS(subhit, file = ff)
  gc()
}

