setwd("~/Desktop/CRUK/britroc-cnsignatures-bfb69cd72c50/manuscript_Rmarkdown")
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