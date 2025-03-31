library(ggbio)
library(GenomicRanges)
library(GenomicAlignments)
k <- 2
bamfile <- BamFile(bamfiles[[k]])
ga <- readGAlignments(bamfile,use.names = TRUE)
p <- autoplot(ga)
print(p)

max.hits.per.sample <- sapply(dflist, function (df) {max(0, max(table(droplevels(df$mrnm))))})
hist(max.hits.per.sample)