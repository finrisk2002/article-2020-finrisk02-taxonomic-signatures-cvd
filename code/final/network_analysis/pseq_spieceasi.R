pseq_spieceasi <- function (pseq, lambda.min.ratio, nlambda, rep.num, seed, thresh) {
	pargs1 <- list(rep.num=rep.num, seed=seed, thresh=thresh)

	res1 <- spiec.easi(pseq, method = "glasso", sel.criterion = "bstars", pulsar.select=TRUE,
					   verbose = TRUE, lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, pargs1)

}