setwd("~/Documents/PhD/other_repos/britroc-cnsignatures-bfb69cd72c50/manuscript_Rmarkdown/")
## Get the data from Geoff
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(QDNAseq))
suppressMessages(library(flexmix))
suppressMessages(library(NMF))
suppressMessages(source("../main_functions.R"))
suppressMessages(source("../helper_functions.R"))
num_cores<-16
samp_annotation<-read.table("data/britroc_sample_data.csv",stringsAsFactors = F,sep = ",",header=T)

result<- samp_annotation %>% 
  filter(!samp_annotation$Failed=="Y") %>%
  dplyr::select(Britroc_No,IM.JBLAB_ID,star_rating) %>%
  dplyr::group_by(Britroc_No) %>%
  dplyr::slice(which.max(star_rating))

all_CN<-readRDS("data/britroc_absolute_copynumber.rds")
all_CN<-all_CN[,colnames(all_CN)%in%samp_annotation[!samp_annotation$Failed=="Y","IM.JBLAB_ID"]]
chrlen<-read.table(paste(this_path,"data/hg19.chrom.sizes.txt",sep="/"),sep="\t",stringsAsFactors = F)[1:24,]

## Use Shiraishi's modified functions
# setwd("~/Desktop/CRUK")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("GenomicRanges")
# install.packages("devtools")
# install.packages("ggplot2")
# install.packages("Rcpp")
library(devtools)
# devtools::install_github("friend1ws/pmsignature")
library(pmsignature)
library(Rcpp)

## Load functions
source("get10Mbfeatures.R")
source("getFeatureVectorAndCountData.R")
source("Turbo_functions.R")
source("SegmentsToSignature.R")

## Get the necessary data to find the signatures
MF0<-get10Mbfeatures(all_CN,chrlen,0)
MFM<-MF0[,c(1,3:9)]

# POFE<-c()
# for(j in 2:8){
#   POFE<-c(POFE,max(MFM[,j])+1)
# }
# FVCD<-getFeatureVectorAndCountData(MFM)
# FEVE<-FVCD[[1]]
# CODA<-FVCD[[2]]
# 
# ## Signature extraction
# SegmentsToSignature(MFM,POFE,FEVE,CODA,3)

##' ---------------------------------------------------------------------
##' ---- Lena's version of modifying the input object of pmsignature with
##' our own segment information

## Read minimal example from pmsignature
inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt.gz", package="pmsignature")
print(inputFile)
G_original <- readMPFile(inputFile, numBases = 5)
# Param_original <- getPMSignature(mutationFeatureData = G_original, K = 3, BG = NULL)

##' we are looking at 5 features. We have 183916 mutations (each with one of these
##' features) across 21 samples. Each mutation has 8 pieces of information, from
##' which the 5 features are derived. We are looking at 3 mutational signatures (k=3).
##' We have these number of compartments in each category: 6 4 4 4 4, which gives 1536 possible
##' combinations. However, only 1534 are observed.


## In our case, our matrix of observations is as follows:
head(MFM)

## we use fewer features for now; I remove the ones with many categories
MFM_small <- MFM[,c(2,7,8)]
head(MFM_small)

## compute the combinations of features
apply(MFM_small, 2, table)

## this can only be done with this function with few (combinations of) categories
# MFM_table <- table(MFM_small)
MFM_small_combinations <- apply(MFM_small, 1, paste0, collapse='-')
num_combinations <- table(MFM_small_combinations)

## modifying the pmsignature object
G_mod <- G_original
names(attributes(G_mod))

##' - featureVectorList
##'   it's a matrix with as many rows as features=5, and entries 1,...,5
##'   it says what are the combinations of features that we can find in our data
##'   sapply(apply(slot(G, "featureVectorList"), 1, table), length)
##'   [1] 6 4 4 4 4
dim(G_mod@featureVectorList)

featureVectorList_mod <- sapply(names(num_combinations), function(j){
  as.numeric(strsplit(j, split = "-")[[1]])})
G_mod@featureVectorList <- featureVectorList_mod
dim(featureVectorList_mod)

##' - countData: it contains the number of mutations in each sample classified in each
##'   category (i.e. combination of features))
##'   it's a matrix with as many rows as observations. The first column contains 1...1534,
##'   i.e. the combinations of features
##'   max(t(slot(G, "countData"))[,1])
##'   [1] 1534
##'   the second column is for the sample id
##'   max(t(slot(G, "countData"))[,2])
##'   [1] 21
##'   the third and last column is the number of mutations in each category
##'   sum(slot(G, "countData")[3,])
##'   [1] 183916
dim(G_mod@countData)
countData_mod0 <- apply(MFM_small, 1, function(j) which(names(num_combinations) == paste0(j, collapse = '-')))
## adding the numerical id of each patient
countData_mod0 <- cbind(patientid=countData_mod0,
                       cat_combination=as.numeric(as.factor(MFM[,1])),
                       count=1)
##' counting the number of instances in each category of combination of features
##' for each patient
countData_mod <- aggregate(count~patientid+cat_combination, countData_mod0, sum)
G_mod@countData <- t(countData_mod)
dim(G_mod@countData)

slotNames(G_mod)
dim(slot(G_mod, "mutationPosition"))
# [1] 183916      8
head(slot(G_mod, "mutationPosition"))
mutationPosition_mod <- cbind.data.frame(chr=NA, pos=NA, ref=NA, alt=NA, strand=NA, context=NA,
                 sampleID=countData_mod0[,'patientid'],
                 mutID=countData_mod0[,'cat_combination'])
slot(G_mod, "mutationPosition") <- mutationPosition_mod

length(slot(G_mod, "sampleList"))
# [1] 21
slot(G_mod, "sampleList") <- unique(MFM$Col0)

(slot(G_mod, "possibleFeatures"))
# [1] 6 4 4 4 4
slot(G_mod, "possibleFeatures") <- apply(MFM_small, 2, function(j) length(unique(j)))

G_mod@type
G_mod@flankingBasesNum
G_mod@transcriptionDirection
G_mod@possibleFeatures

max(slot(G_mod, "mutationPosition")$sampleID)

# save.image(file = "~/Desktop/image_Jorge.RData")
load("~/Desktop/image_Jorge.RData")

# results_signatures <- getPMSignature(mutationFeatureData = G_mod, K = 3, BG = NULL)
# results_signatures@signatureFeatureDistribution

# image(results_signatures@sampleSignatureDistribution)
# 
# 
# dim(G_original@featureVectorList)
# length(G_original@sampleList)
# dim(G_original@countData)
# dim(G_original@mutationPosition)
# (G_original@type)
# (G_original@flankingBasesNum)
# (G_original@transcriptionDirection)
# (G_original@possibleFeatures)
# 
# dim(G_mod@featureVectorList)
# length(G_mod@sampleList)
# dim(G_mod@countData)
# dim(G_mod@mutationPosition)
# (G_mod@type)
# (G_mod@flankingBasesNum)
# (G_mod@transcriptionDirection)
# (G_mod@possibleFeatures)
# 
# 
# prod(G_original@possibleFeatures)
# prod(G_mod@possibleFeatures)


##' getPMSignature(mutationFeatureData = G_mod, K = 3, BG = NULL) makes R crash so
##' here I am trying to debug it

mutationFeatureData = G_mod
K = 3
BG = NULL
numInit = 10
tol = 1e-04
maxIter = 10000
# function (mutationFeatureData, K, BG = NULL, numInit = 10, tol = 1e-04, 
#           maxIter = 10000) 
# {
  if (!is.null(BG)) {
    isBG <- TRUE
    varK <- K - 1
  }else {
    isBG <- FALSE
    varK <- K
    BG <- 0
  }
  sampleNum <- length(slot(mutationFeatureData, "sampleList"))
  fdim <- slot(mutationFeatureData, "possibleFeatures")
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
    p0 <- c(pmsignature:::convertToTurbo_F(as.vector(F), fdim, K, isBG), 
            pmsignature:::convertToTurbo_Q(as.vector(t(Q)), K, sampleNum))
    Y <- list(list(sampleNum, fdim, slot(mutationFeatureData, 
                                         "featureVectorList"), slot(mutationFeatureData, "countData")), 
              K, isBG, BG)
    res1 <- pmsignature:::mySquareEM(p0, Y, tol = tol, maxIter = maxIter)
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
  return(new(Class = "EstimatedParameters", type = slot(mutationFeatureData, 
                                                        "type"), flankingBasesNum = slot(mutationFeatureData, 
                                                                                         "flankingBasesNum"), transcriptionDirection = slot(mutationFeatureData, 
                                                                                                                                            "transcriptionDirection"), possibleFeatures = slot(mutationFeatureData, 
                                                                                                                                                                                               "possibleFeatures"), sampleList = slot(mutationFeatureData, 
                                                                                                                                                                                                                                      "sampleList"), signatureNum = as.integer(K), isBackGround = isBG, 
             backGroundProb = BG, signatureFeatureDistribution = F, 
             sampleSignatureDistribution = Q, loglikelihood = tempL))
}
<bytecode: 0x7ffe3e289350>
  <environment: namespace:pmsignature>
  
  
  pmsignature:::mySquareEM(p0, Y, tol = tol, maxIter = maxIter)
p=p0
y=Y
tol = tol
maxIter = maxIter
# function (p, y, tol = 1e-04, maxIter = 10000) 
# {
  prevL <- -Inf
  step.min <- 1
  step.max <- 1
  step.max0 <- 1
  mstep <- 4
  objfn.inc <- 1
  updEvalNum <- 0
  LEvalNum <- 0
  useSquareEM <- 0
  iterNum <- 0
  convFlag <- FALSE
  startTime <- proc.time()
  newL <- pmsignature:::calcPMSLikelihood(p, y)
  LEvalNum <- LEvalNum + 1
  for (iterNum in 1:maxIter) {
    p1 <- pmsignature::updatePMSParam(p, y)
    updEvalNum <- updEvalNum + 1
    if (any(is.nan(unlist(p1)))) {
      stop("Error in function evaluation")
    }
    q1 <- p1 - p
    sr2 <- crossprod(q1)
    p2 <- updatePMSParam(p1, y)
    updEvalNum <- updEvalNum + 1
    if (any(is.nan(unlist(p2)))) {
      stop("Error in function evaluation")
    }
    q2 <- p2 - p1
    sq2 <- sqrt(crossprod(q2))
    sv2 <- crossprod(q2 - q1)
    srv <- crossprod(q1, q2 - q1)
    alpha <- -srv/sv2
    alpha <- max(step.min, min(step.max, alpha))
    p.new <- p + 2 * alpha * q1 + alpha^2 * (q2 - q1)
    if (isTRUE(abs(alpha - 1) > 0.01)) {
      p.new <- updatePMSParam(p.new, y)
      updEvalNum <- updEvalNum + 1
    }
    if (any(is.nan(p.new)) | !PMSboundary(y)(p.new)) {
      p.new <- p2
      newL <- calcPMSLikelihood(p2, y)
      LEvalNum <- LEvalNum + 1
      if (isTRUE(all.equal(alpha, step.max))) {
        step.max <- max(step.max0, step.max/mstep)
      }
      alpha <- 1
    }
    else {
      newL <- calcPMSLikelihood(p.new, y)
      LEvalNum <- LEvalNum + 1
      if (is.nan(newL) | (newL > prevL + objfn.inc)) {
        p.new <- p2
        lnew <- calcPMSLikelihood(p2, y)
        LEvalNum <- LEvalNum + 1
        if (alpha == step.max) {
          step.max <- max(step.max0, step.max/mstep)
        }
        alpha <- 1
      }
      else {
        useSquareEM <- useSquareEM + 1
      }
    }
    if (isTRUE(all.equal(alpha, step.max))) {
      step.max <- mstep * step.max
    }
    if (step.min < 0 & isTRUE(all.equal(alpha, step.min))) {
      step.min <- mstep * step.min
    }
    p <- p.new
    if (abs(prevL - newL) < tol) {
      convFlag <- TRUE
      break
    }
    if (!is.nan(newL)) {
      prevL <- newL
    }
  }
  calcTime <- proc.time() - startTime
  return(list(par = p, value.objfn = newL, itr = iterNum, fpeval = updEvalNum, 
              convergence = convFlag, elapsed.time = calcTime[3]))
}