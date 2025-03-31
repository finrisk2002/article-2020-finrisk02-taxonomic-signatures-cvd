library(biomformat)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(viridis)


plot_pcoa <- function(ord.pcoa, ord.df, titlestr) {
  p <- ggplot(data=ord.df, aes(x=pcoa1, y=pcoa2)) + geom_point(size=0.4, alpha=0.4)
  p <- p + ggtitle(titlestr)
  eig.fracs <- ord.pcoa$values$Relative_eig
  p <- p + xlab(paste("PCoA axis 1,", round(100*eig.fracs[1],1), "%"))
  p <- p + ylab(paste("PCoA axis 2,", round(100*eig.fracs[2],1), "%"))
  p
}

plot_pcoa_abudcolor <- function (ord.pcoa, ord.df, titlestr, colorvar) {

  p <- ggplot(data=ord.df, aes_string(x="pcoa1", y="pcoa2", colour=colorvar)) + geom_point(size=0.4)
  p <- p+ scale_color_distiller(palette='Spectral', name='Rel. abundance')
  p <- p + ggtitle(titlestr)
  eig.fracs <- ord.pcoa$values$Relative_eig
  p <- p + xlab(paste("PCoA axis 1,", round(100*eig.fracs[1],1), "%"))
  p <- p + ylab(paste("PCoA axis 2,", round(100*eig.fracs[2],1), "%"))
}
