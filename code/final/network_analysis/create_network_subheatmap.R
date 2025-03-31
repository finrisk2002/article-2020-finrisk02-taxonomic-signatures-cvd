# heatmap of subnetwork "activations"
create_network_subheatmap <- function (n.tmp.df, order.samples, order.components,
  n.components, value="value", n=NULL, text.size=15, legend.title.size=15) {

  if (!is.null(n)) {
    n.tmp.df <- filter(n.tmp.df, OTU %in% (n %v% "taxa.names"))
  }
  n.tmp.df$OTU <- as.character(n.tmp.df$OTU)

  n.tmp.df[["SampleID"]] <- factor(n.tmp.df[["SampleID"]], levels=order.samples)

  if (!("ClusterNo" %in% colnames(n.tmp.df))) {
    n.tmp.df[['ClusterNo']] <- sapply(n.tmp.df[['OTU']], function (otu) {n.components[otu,"componentID"]})
  }

  n.tmp.df[["ClusterNo"]] <- factor(n.tmp.df[["ClusterNo"]], ordered=TRUE)
  #n.tmp.df[["ClusterName"]] <- sapply(n.tmp.df[["ClusterNo"]],
  #                                                function (x) {get_cluster_name(n.components, x)})

  n.tmp.df[["OTU.Subnet"]] <- sapply(1:nrow(n.tmp.df), function(ri) {paste(n.tmp.df[ri,"OTU"], n.tmp.df[ri,"ClusterNo"])})
  n.tmp.df[["OTU"]] <- factor(n.tmp.df[["OTU"]], levels=order.components)
  n.tmp.df[["OTU.Subnet"]] <- factor(n.tmp.df[["OTU.Subnet"]], levels=order.components)

  n.tmp.df[["value.clipped"]] <- n.tmp.df[[value]]
  n.tmp.df[n.tmp.df[["value.clipped"]] > 5,"value.clipped"] <- 5
  n.tmp.df[n.tmp.df[["value.clipped"]] < -5,"value.clipped"] <- -5

  p <- ggplot(n.tmp.df, aes_string(x="SampleID", y="OTUnice", fill="value.clipped")) + geom_tile()

  tmpabsmax <- signif(1.05*max(abs(n.tmp.df[[value]])), digits=2)
  tmpmax <- signif(1.05*max(n.tmp.df[[value]]), digits=2)
  tmpmin <- signif(1.05*min(n.tmp.df[[value]]), digits=2)


  p <- p + scale_fill_gradientn(name="Abundance",
                                colours=c("darkblue", "blue","white","red","darkred"),
                                values=rescale(c(-5,-2,0,2,5)), limits=c(-5,5), breaks=c(-5,-2.5, 0, 2.5, 5),
                                labels=c("< -5", "-2.5", "0", "2.5", "5 <"))
  p <- p + xlab("") + ylab("")
  p <- p + theme_classic(text.size)
  p <- p + theme(axis.text.x = element_blank())
  p <- p + theme(axis.text.y = element_text(hjust = 1.0, size=floor(1.0*text.size), color="black"))
  p <- p + theme(axis.title.x = element_blank())
  p <- p + theme(axis.title.y = element_blank())

  p <- p + theme(axis.ticks.x = element_blank())
  p <- p + theme(axis.line.y = element_blank())
  p <- p + theme(axis.line.x = element_blank())
  p <- p + theme(plot.margin = unit(c(0.1,0.1,0.1,0.5), "cm"))
  p <- p + theme(legend.title = element_text(size=legend.title.size))
  p <- p + theme(legend.text = element_text(size=text.size))
  p

}