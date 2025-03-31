library(Rsamtools)
source("funcs.R")

#which <- IRangesList(seq1=IRanges(1000, 2000),
#                     seq2=IRanges(c(100, 1000),#
#		     c(1000, 2000)))
#param <- ScanBamParam(which=which, what=what)

#what <- c("qname","flag", "rname","strand", "pos", "qwidth", "mapq", "cigar","mrnm", "mpos", "isize","seq",  "qual")

what <- c("qname","flag", "rname","strand", "pos", "qwidth", "mapq", "cigar","mpos", "isize")

param <- ScanBamParam(what=what)

bamfiles <- list.files("~/data/tmp/BAM",
                      pattern = ".bam$",
                      full.names = TRUE)


# Only keep the samples that are in FR02
bamfiles <- bamfiles[gsub(".bam$", "", gsub("/home/lemila/data/tmp/BAM/", "", bamfiles)) %in% include.samples]

dflist <- list()
# bamfiles <- setdiff(bamfiles, discard)
vfdb <- readLines("~/data/FINRISK/VFDB/VFDB_setB_nt.fas")
headers <- vfdb[!grepl("^[A-Z]", vfdb)]
headers <- headers[sapply(headers, nchar) > 0]
headers <- gsub("^>", "", headers)

for (bamfile in bamfiles) {

  print(paste(bamfile, match(bamfile, bamfiles), "/", length(bamfiles)))
    
  bam <- scanBam(bamfile, param=param)
  #bam <- unname(bam) # names not useful in unlisted result  
  elts <- setNames(bamWhat(param), bamWhat(param))
  lst <- lapply(elts, function(elt) .unlist(lapply(bam, "[[", elt)))
  df <- do.call("DataFrame", lst)
  # df$bam <- rep(bamfile, nrow(df))

  # Retrieve the full header based on rname
  strs <- as.character(df$rname)
  strs.uniq <- unique(strs)
  h <- unname(sapply(strs.uniq, function (s) unique(headers[grepl(s, headers, fixed = TRUE)])))
  h <- h[match(strs, strs.uniq)]
  if (!length(h) == nrow(df)) {stop("check")}
  df$rname_header <- h

  dflist[[bamfile]] <- df

}

print("BAM files read ok")

# Just keep the IDs
library(stringr)
names(dflist) <- word(word(names(dflist), sep = "/", 7), sep = "\\.", 1)
saveRDS(dflist, file = "dflist.Rds")

# dfl <- do.call("rbind", dflist)


