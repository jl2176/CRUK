## DESCRIPTION

#This is a modified version of the function getPMSignature from the package
#pmsignature. Inputs:
#Segments_matrix <- Corresponds to the output from "get10Mbfeatures"
#Possible_features <- Vector with the number of possible values for each feature
#Feature_vector <- Matrix of dimensions: (number of features) x (number of observed combinations of features) containing the observed combinations of features
#Count_data <- Matrix of dimensions 3 rows, each containing: id of combination of features, sample id and number of repetitions. 
#Note that if a combination of features is not found in a sample the column is not included
#K <- Number of mutational signatures


SegmentsToSignature<-function (Segments_matrix, Possible_features, Feature_Vector, Count_Data, K, BG = NULL, numInit = 10, tol = 1e-04, 
          maxIter = 10000) 
{
  if (!is.null(BG)) {
    isBG <- TRUE
    varK <- K - 1
  }
  else {
    isBG <- FALSE
    varK <- K
    BG <- 0
  }
  sampleNum <- length(unique(Segments_matrix[,1]))
  fdim <- Possible_features
  tempL <- -Inf
  tempPar <- c()
  for (kkk in 1:numInit) {
    F <- array(0, c(varK, length(fdim), max(fdim)))
    for (k in 1:varK) {
      for (kk in 1:length(fdim)) {
        F[k, kk, 1:fdim[kk]] <- rgamma(fdim[kk], rep(1, 
                                                     fdim[kk]))
        F[k, kk, 1:fdim[kk]] <- F[k, kk, 1:fdim[kk]]/sum(F[k, 
                                                           kk, 1:fdim[kk]])
      }
    }
    Q <- matrix(rgamma(sampleNum * K, 1, 1), K, sampleNum)
    Q <- sweep(Q, 2, apply(Q, 2, sum), `/`)
    p0 <- c(convertToTurbo_F(as.vector(F), fdim, K, isBG), 
            convertToTurbo_Q(as.vector(t(Q)), K, sampleNum))
    Y <- list(list(sampleNum, fdim, Feature_Vector, Count_Data), 
              K, isBG, BG)
    res1 <- mySquareEM(p0, Y, tol = tol, maxIter = maxIter)
    cat(paste("#trial: ", sprintf("%2d", kkk), "; #iteration: ", 
              sprintf("%4d", as.integer(res1$itr)), "; time(s): ", 
              sprintf("%4.2f", res1$elapsed.time), "; convergence: ", 
              res1$convergence, "; loglikelihood: ", sprintf("%.4f", 
                                                             res1$value.objfn), "\n", sep = ""))
    if (res1$value.objfn > tempL) {
      tempL <- res1$value.objfn
      tempPar <- res1$par
    }
  }
  lenF <- varK * (sum(fdim) - length(fdim))
  lenQ <- sampleNum * (K - 1)
  F <- convertFromTurbo_F(tempPar[1:lenF], fdim, K, isBG)
  Q <- convertFromTurbo_Q(tempPar[(lenF + 1):(lenF + lenQ)], 
                          K, sampleNum)
  dim(F) <- c(varK, length(fdim), max(fdim))
  dim(Q) <- c(sampleNum, K)
  #type, flanking bases and transcription direction
  return(new(Class = "EstimatedParameters", possibleFeatures = Possible_features, sampleList = unique(Segments_matrix[,1]), signatureNum = as.integer(K), isBackGround = isBG, backGroundProb = BG, signatureFeatureDistribution = F, sampleSignatureDistribution = Q, loglikelihood = tempL))
}