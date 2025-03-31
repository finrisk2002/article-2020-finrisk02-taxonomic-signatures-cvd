library(stringr)
library(glue)
#path <- "~/tmp/data"
path <- "data"
fs <- list.files(path, full.names = TRUE)
#ids <- unique(word(word(fs, sep = "\\."), sep = "/", start = 2))
ids <- unique(str_replace_all(fs, "data/|.R[1|2].trimmed.filtered.fastq.gz", ""))

bamfiles <- list.files("~/data/tmp/BAM/", pattern = ".bam$")
ids <- setdiff(ids, gsub(".bam", "", bamfiles))
ids <- setdiff(ids, "files.txt")

conn <- file("samples.txt", "w")
writeLines(paste("sample", "r1", "r2", sep = "\t"), con = conn)
for (id in ids) {
  writeLines(paste(glue("{id}"),
                   glue("{id}.R1.trimmed.filtered.fastq.gz"),
		   glue("{id}.R2.trimmed.filtered.fastq.gz"),
		   sep = "\t"), con = conn)
}

close(conn)
