# SpiecEasi results to edgelist format, where edgelist will containt useful variables for plots
get_spiec_edgelist <- function(spiec.res, spiec.res.a, pseq.core, subnetthr,
                                            edgeWcoef=20, edgeWmin=1) {

  tmp.edglist <- filter(melt(spiec.res.a), value > 0)
  tmp.edglist$Var1 <- as.character(tmp.edglist$Var1)
  tmp.edglist$Var2 <- as.character(tmp.edglist$Var2)

  tmp.edglist <- tmp.edglist %>% transmute(Var1=Var1, Var2=Var2, SpiecEasiNotZero = value)

  # correlation matrix, also to edgelist
  spiec.res.c  <- cov2cor(as.matrix(getOptCov(spiec.res)))
  diag(spiec.res.c) <- 0
  spiec.res.c[lower.tri(spiec.res.c)] <- 0
  colnames(spiec.res.c) <- rownames(spiec.res.c) <- colnames(t(abundances(pseq.core)))
  tmp <- melt(spiec.res.c)
  tmp$Var1 <- as.character(tmp$Var1)
  tmp$Var2 <- as.character(tmp$Var2)
  tmp <- tmp %>% transmute(Var1=Var1, Var2=Var2, SpiecEasiEdge = value)

  # merge correlations to main edgelist
  tmp.edglist <- merge(tmp.edglist, tmp, by=c("Var1", "Var2"))
  tmp.edglist$SpiecEasiEdgeAbs <- abs(tmp.edglist$SpiecEasiEdge)
  tmp.edglist$SpiecEasiEdgePr <- sapply(round(tmp.edglist$SpiecEasiEdge, digits=5),
                                        function(x) {sprintf("%0.5f",x)})
  tmp.edglist$SpiecEasiEdgeW <- pmax(edgeWmin, edgeWcoef*abs(tmp.edglist$SpiecEasiEdge))

  tmp.edglist <- filter(tmp.edglist, SpiecEasiEdgeAbs>subnetthr)


  tmp.edglist
}