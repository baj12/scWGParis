# filter a sparse expression matrix
# genes in rows x cells in column

filterSparseMat <- function(sparseMat, minCellPerGene=3, minGenePerCell=1000, maxGenePerCell=10000){

	filtMat <- sparseMat <- sparseMat[ rowSums(sparseMat) > minCellPerGene, ]
	filtMat <- filtMat[ , colSums(filtMat) > minGenePerCell & colSums(filtMat) < maxGenePerCell ]

	return(filtMat)

}


# ensuire a list of matrices contains the same genes
# required before mnn correction

adjustMatContent <- function(listOfMatrices){

	sharedGenes <- Reduce(intersect, lapply(listOfMatrices, rownames))
	listOfEqMat <- lapply(listOfMatrices, function(x) {x[sharedGenes,]})

	return(listOfEqMat)

}


# evaluate the batch effect with kBET
# see
# https://github.com/theislab/kBET
# it unveals difference only for batch3

testBatch_kBET <- function(matrice1, matrice2){

	require(kBET)
	b_est <- kBET((as.matrix(cbind(matrice1, matrice2))), c(rep('batch1', ncol(matrice1)), rep('batch2', ncol(matrice2))), plot=FALSE)
	
	return(b_est$average.pval)

}
