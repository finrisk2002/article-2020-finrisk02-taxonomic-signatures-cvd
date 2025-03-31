# hclust on correlation distance
get_cor_hclust_order <- function(transformed_abundance_matrix) {
  mat<- transformed_abundance_matrix
  scor.mat <- cor(mat, use="pairwise.complete.obs")
  scor.mat[is.na(scor.mat)] <-0
  rind <- hclust(as.dist(1 - scor.mat))$order
  order.items <- colnames(mat)[rind]

  order.items
}