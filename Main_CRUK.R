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
setwd("~/Desktop/CRUK")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
install.packages("devtools")
install.packages("ggplot2")
install.packages("Rcpp")
library(devtools)
devtools::install_github("friend1ws/pmsignature")
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
POFE<-c()
for(j in 2:8){
  POFE<-c(POFE,max(MFM[,j])+1)
}
FVCD<-getFeatureVectorAndCountData(MFM)
FEVE<-FVCD[[1]]
CODA<-FVCD[[2]]

## Signature extraction
SegmentsToSignature(MFM,POFE,FEVE,CODA,3)

