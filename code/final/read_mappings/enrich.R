library(dplyr)
library(microbiome)

contig.length.threshold <- 500
threshold.coverage <- 0.9
accepted.hits <- stat %>%
         filter(max.contig.length >= contig.length.threshold) %>% # At least one continuous 500bp hit
         filter(coverage.rel >= threshold.coverage) # High coverage for the gene	 



esc <- abundances(microbiome::transform(gen, "clr"))["g_Escherichia", ]

enrichment.pvalues <- list()

# Are high-Escherichia samples enriched with specific virulence genes?
pvs <- c()
M <- meta(gen)
for (g in unique(accepted.hits[, "gene"])) {
  print(g)
  # g <- "VFG049154(gb|YP_006635480.1)"          # Selected VFDB gene
  s <- unique(subset(accepted.hits, gene == g)$sample) # Samples matching this gene
  # sc <- setdiff(names(esc), s) # Complement
  # Grouping: samples that match and do not match this gene
  gr <- factor(names(esc) %in% s)  
  # Enrichment p-value by logistic regression (Eschericia abundances predicting gene presence)

  pv <- coef(summary(glm(gr ~ esc + M$BL_AGE + M$BMI + M$MEN + M$CURR_SMOKE + M$PREVAL_DIAB + M$SYSTM + M$BL_USE_RX_L + M$BP_TREAT, family = binomial(link = "logit"))))["esc", "Pr(>|z|)"] 
  # boxplot(esc ~ gr, main = paste(g, pv))
  pvs[[g]] <- pv
}
pva <- p.adjust(pvs, "BH")
enrichment.pvalues[["Escherichia"]] <- pva


# Living vs. Dead people: differences in virulence genes?
pvs2 <- c()
for (g in unique(accepted.hits[, "gene"])) {
  print(g)
  s <- unique(subset(accepted.hits, gene == g)$sample) # Samples matching this gene
  # Grouping: samples that match and do not match this gene
  gr <- rownames(M) %in% s + 0
  # Enrichment p-value by logistic regression (Eschericia abundances predicting gene presence)
  M$group <- as.numeric(gr)
  D <- M[, c("group", "DEATH", "BL_AGE", "BMI", "MEN", "CURR_SMOKE", "PREVAL_DIAB", "SYSTM", "BL_USE_RX_L", "BP_TREAT")]

  D <- D[rowSums(is.na(D)) == 0, ]
  pv <- coef(summary(glm(group ~ DEATH + BL_AGE + BMI + MEN + CURR_SMOKE + PREVAL_DIAB + SYSTM + BL_USE_RX_L + BP_TREAT, data = D, family = binomial(link = "logit"))))["DEATH", "Pr(>|z|)"] 
  pvs2[[g]] <- pv
}
pva2 <- p.adjust(pvs2, "BH")
enrichment.pvalues[["Death"]] <- pva2




