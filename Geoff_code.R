## Load packages and other functions
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(QDNAseq))
suppressMessages(library(flexmix))
suppressMessages(library(NMF))
suppressMessages(source("../main_functions.R"))
suppressMessages(source("../helper_functions.R"))
num_cores<-16

## Set ggplot theme
cbPalette <- c(RColorBrewer::brewer.pal(8,"Dark2"),RColorBrewer::brewer.pal(9,"Set1"),"black")
ggplot2::theme_set(ggplot2::theme_gray(base_size = 5))
my_theme<-ggplot2::theme_bw()+ggplot2::theme(axis.text=ggplot2::element_text(size=5),
                                             axis.title=ggplot2::element_text(size=5),
                                             strip.text.x = ggplot2::element_text(size = 7),
                                             strip.text.y = ggplot2::element_text(size = 7),
                                             legend.text = ggplot2::element_text(size = 7),
                                             panel.grid.minor = ggplot2::element_blank(),
                                             panel.grid.major = ggplot2::element_blank())
#Sample annotation
samp_annotation<-read.table("data/britroc_sample_data.csv",stringsAsFactors = F,sep = ",",header=T)

table(samp_annotation[samp_annotation$Failed=="Y","Notes"])

table(samp_annotation[!samp_annotation$Failed=="Y","star_rating"])

result<- samp_annotation %>% 
  filter(!samp_annotation$Failed=="Y") %>%
  dplyr::select(Britroc_No,IM.JBLAB_ID,star_rating) %>%
  dplyr::group_by(Britroc_No) %>%
  dplyr::slice(which.max(star_rating))

## WGS ploidy purity correlation
# Read in ploidy and purity estimates obtained from Battenberg plots.

deep_metrics <- read.csv("data/britroc_deepWGS_ploidy_purity_49.csv", header=T, as.is = T)
deep_metrics$purity <- deep_metrics$purity/100
deep_metrics <- deep_metrics %>%
  mutate(name= sub("[0-9]+_tumor_","",x=wgs_ID)) %>%
  dplyr::select(name,ploidy, purity) %>%
  tidyr::gather("metric", "deep", -name)

# Get ploidy for sWGS samples

CN <- readRDS("data/britroc_absolute_copynumber.rds")
Biobase::pData(CN)$ploidy <- getPloidy(CN) %>%
  .$out

shallow_metrics <- Biobase::pData(CN) %>%
  dplyr::select(name, ploidy, purity) %>%
  filter(name %in% deep_metrics$name) %>%
  tidyr::gather("metric", "shallow", -name)


# Annotate sWGS samples with star_rating
samp_annotation_all <- read.csv("data/britroc_sample_data.csv", as.is=T)
samp_annot <- samp_annotation_all %>% 
  filter(IM.JBLAB_ID %in% shallow_metrics$name & Failed !="Y") %>%
  dplyr::select(IM.JBLAB_ID , star_rating)

# Combine estimates from deep and shallow WGS methods and retain 3-star samples
metrics <- inner_join(deep_metrics, shallow_metrics, by = c("name", "metric")) %>%
  inner_join(samp_annot,c("name" = "IM.JBLAB_ID" ) ) %>%
  filter(star_rating ==3)

# Perform correlation test of ploidy, purity in Britroc samples with dWGS and sWGS data
correlations <- metrics %>% 
  group_by(metric) %>%
  do(broom::tidy(cor.test(.$deep,.$shallow)))
# Create dataframes for annotation and generate plot for 3-star samples.

correlations$metric <- factor(correlations$metric, levels=c("ploidy", "purity"))

annot_text <- correlations %>% 
  ungroup() %>%
  dplyr::select(metric, estimate, p.value) %>%
  mutate(x=c(2,0.4),y.cor=c(4.5,0.9), y.pval= c(4.375, 0.87)) %>%
  group_by(metric) %>%
  mutate(label.cor=paste0("R^2== ",round(estimate,2)), label.pval = paste0("P== ", format(signif(p.value,1),digits=1, scientific = T)) )


cor_plot <- metrics %>%
  ggplot(aes(x=deep, y=shallow)) +
  geom_point()  + 
  geom_smooth(method = "lm") + 
  facet_wrap(~metric, scales = "free") +
  labs(x= "Deep WGS", y= "sWGS") +
  theme_bw()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=7),
        strip.text.x = element_text(size = 7), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


cor_plot + 
  geom_text(data=annot_text, aes(label=label.cor, x=x,y=y.cor), parse = T,inherit.aes=FALSE, size=2) +
  geom_text(data=annot_text, aes(label=label.pval, x=x,y=y.pval), parse = T,inherit.aes=FALSE, size=2)

## Back to Signature identification validation

#read in absolute CN post ABCEL processing
all_CN<-readRDS("data/britroc_absolute_copynumber.rds")
all_CN<-all_CN[,colnames(all_CN)%in%samp_annotation[!samp_annotation$Failed=="Y","IM.JBLAB_ID"]]

#extract 3 star CN from each case
ids<-result %>% filter(star_rating==3)
ids<-ids$IM.JBLAB_ID
hq_CN<-all_CN[,colnames(all_CN)%in%ids]
ids<-result %>% filter(star_rating==2)
ids<-ids$IM.JBLAB_ID
lq_CN<-all_CN[,colnames(all_CN)%in%ids]

#estract copy-number features
CN_features<-extractCopynumberFeatures(hq_CN,cores=num_cores)


## Apply mixture models
seed=77777
min_prior=0.001
model_selection="BIC"
nrep=1
niter=1000

dat<-as.numeric(CN_features[["segsize"]][,2])
segsize_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                         min_prior=min_prior,niter=niter,nrep=nrep,min_comp=10,max_comp=10)

dat<-as.numeric(CN_features[["bp10MB"]][,2])
bp10MB_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                        min_prior=min_prior,niter=niter,nrep=nrep,min_comp=3,max_comp=3)

dat<-as.numeric(CN_features[["osCN"]][,2])
osCN_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=3,max_comp=3)

dat<-as.numeric(CN_features[["bpchrarm"]][,2])
bpchrarm_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                          min_prior=min_prior,niter=niter,nrep,min_comp=2,max_comp=5)

dat<-as.numeric(CN_features[["changepoint"]][,2])
changepoint_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                             min_prior=min_prior,niter=niter,nrep=nrep,min_comp=7,max_comp=7)

dat<-as.numeric(CN_features[["copynumber"]][,2])
copynumber_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                            nrep=nrep,min_comp=2,max_comp=10,min_prior=0.005,niter=2000)

CN_components<-list(segsize=segsize_mm,bp10MB=bp10MB_mm,osCN=osCN_mm,changepoint=changepoint_mm,copynumber=copynumber_mm,bpchrarm=bpchrarm_mm)

britroc_sample_component_matrix<-generateSampleByComponentMatrix(CN_features,CN_components,cores=1,subcores=num_cores)
NMF::aheatmap(britroc_sample_component_matrix,fontsize = 7,Rowv=FALSE,Colv=FALSE,legend = T,breaks=c(seq(0,199,2),500),main="Component x Sample matrix")

## Signature selection (THIS SELECTS THE OPTIMAL NUMBER OF SIGNATURES)
nmfalg<-"brunet"
# seed<-77777
# estim.r <- NMF::nmfEstimateRank(t(britroc_sample_component_matrix), 3:12,seed = seed,nrun=1000,
#                                 verbose=F,method=nmfalg,.opt = paste0("p",num_cores))
# V.random <- randomize(t(britroc_sample_component_matrix))
# estim.r.random <- NMF::nmfEstimateRank(V.random, 3:12, seed =seed,nrun=1000,
#                                        verbose=F,method=nmfalg,.opt = paste0("p",num_cores))
# p<-plot(estim.r,estim.r.random, 
#         what = c("cophenetic", "dispersion","sparseness", "silhouette"),xname="Observed",yname="Randomised",main="")+
#   theme(axis.text=element_text(size=5),axis.title=element_text(size=5),
#         strip.text.x = element_text(size = 5),
#         strip.text.y = element_text(size = 5),
#         legend.text = element_text(size = 5),
#         legend.title = element_text(size = 7))
# g<-ggplotGrob(p)
# g[["grobs"]][[2]]$children[[4]]$size[[1]]<-0.5
# g[["grobs"]][[3]]$children[[4]]$size[[1]]<-0.5
# g[["grobs"]][[4]]$children[[4]]$size[[1]]<-0.5
# g[["grobs"]][[5]]$children[[4]]$size[[1]]<-0.5
# grid::grid.newpage()
# grid::grid.draw(g)

nsig<-7
sigs<-NMF::nmf(t(britroc_sample_component_matrix),nsig,seed=seed,nrun=1000,method=nmfalg,.opt = paste0("p",num_cores))
png('PatientSignature5.png')
pdf(file = 'PxS.pdf')
coefmap(sigs,Colv="consensus",tracks=c("basis:"), main="Patient x Signature matrix")
dev.off()

basismap(sigs,Rowv=NA,main="Signature x Component matrix")


### PCAWG
pcawg_CN_features<-readRDS("data/pcawg_CN_features.rds")
pcawg_sample_component_matrix<-generateSampleByComponentMatrix(pcawg_CN_features,CN_components)

NMF::aheatmap(pcawg_sample_component_matrix,Rowv=NULL, main="Component x Sample matrix")


## Survival analysis
#britroc feat_sig matrix
feat_sig_mat<-basis(sigs)
reord_britroc<-as.integer(c(2,6,5,4,7,3,1))
names(reord_britroc)<-paste0("s",1:7)
feat_sig_mat<-feat_sig_mat[,reord_britroc]
colnames(feat_sig_mat)<-paste0("s",1:nsig)
sig_feat_mat<-t(feat_sig_mat)

sig_thresh<-0.01

sig_pat_mat_hq<-scoef(sigs)
sig_pat_mat_hq<-sig_pat_mat_hq[reord_britroc,]
rownames(sig_pat_mat_hq)<-paste0("s",1:nsig)
sig_pat_mat_hq<-normaliseMatrix(sig_pat_mat_hq,sig_thresh)

sig_pat_mat_pcawg<-scoef(pcawg_sigs)
rownames(sig_pat_mat_pcawg)<-paste0("s",1:nsig)
sig_pat_mat_pcawg<-sig_pat_mat_pcawg[reord_pcawg,]
rownames(sig_pat_mat_pcawg)<-paste0("s",1:nsig)
sig_pat_mat_pcawg<-normaliseMatrix(sig_pat_mat_pcawg,sig_thresh)

sig_pat_mat_tcga<-scoef(tcga_sigs)
rownames(sig_pat_mat_tcga)<-paste0("s",1:nsig)
sig_pat_mat_tcga<-sig_pat_mat_tcga[reord_tcga,]
rownames(sig_pat_mat_tcga)<-paste0("s",1:nsig)
sig_pat_mat_tcga<-normaliseMatrix(sig_pat_mat_tcga,sig_thresh)

#extract features for lower quality samples
lowquality_britroc_CN_features<-extractCopynumberFeatures(lq_CN)
lowquality_britroc_sample_component_matrix<-generateSampleByComponentMatrix(lowquality_britroc_CN_features,CN_components)
britroc_lq_ids<-rownames(lowquality_britroc_sample_component_matrix)

#assign signatures
sig_pat_mat_lq<-YAPSA::LCD(t(lowquality_britroc_sample_component_matrix),feat_sig_mat)
rownames(sig_pat_mat_lq)<-paste0("s",1:nsig)
sig_pat_mat_lq<-normaliseMatrix(sig_pat_mat_lq,sig_thresh)

#assign signatures to remaining 2 and 3 star samples (for cases with multiple samples)
remain_samp<-samp_annotation[(!samp_annotation$IM.JBLAB_ID%in%colnames(cbind(sig_pat_mat_hq,sig_pat_mat_lq)))&(!samp_annotation$star_rating==1),]
remain_CN<-all_CN[,colnames(all_CN)%in%remain_samp$IM.JBLAB_ID]

remain_britroc_CN_features<-extractCopynumberFeatures(remain_CN)
remain_britroc_sample_component_matrix<-generateSampleByComponentMatrix(remain_britroc_CN_features,CN_components)

britroc_remain_ids<-rownames(remain_britroc_sample_component_matrix)

sig_pat_mat_remain<-YAPSA::LCD(t(remain_britroc_sample_component_matrix),feat_sig_mat)
rownames(sig_pat_mat_remain)<-paste0("s",1:nsig)
sig_pat_mat_remain<-normaliseMatrix(sig_pat_mat_remain,sig_thresh)

sig_pat_mat_britroc<-cbind(sig_pat_mat_hq,sig_pat_mat_lq)
sig_pat_mat_britroc_all<-cbind(sig_pat_mat_britroc,sig_pat_mat_remain)



library(survival)
set.seed(seed)
samp_num<-0.7 #fraction of the data considered training

tcga_clin<-read.table("data/tcga_sample_info.tsv",stringsAsFactors = F,header=T,sep="\t")
pcawg_clin<-read.table("data/pcawg_sample_info.tsv",sep="\t",header=T,stringsAsFactors = F)
pcawg_cohort<-read.table("data/pcawg_cohort_info.tsv",header=T,sep="\t",stringsAsFactors = F)

#britroc survival data
pat_info<-samp_annotation[,c("Britroc_No","IM.JBLAB_ID","star_rating")]
pat_info<-samp_annotation[grepl("IM",samp_annotation$IM.JBLAB_ID),]
pat_info<-pat_info[order(pat_info$Britroc_No,pat_info$star_rating,decreasing = T),]
pat_info<-pat_info[!duplicated(pat_info$Britroc_No),]
pat_info<-pat_info[!pat_info$IM.JBLAB_ID%in%c("IM_91","IM_70"),]#remove misclassified samples
surv_dat<-read.table("data/britroc_survival_intervals.tsv",sep="\t",header=T,stringsAsFactors = F)
surv_dat<-surv_dat[!duplicated(surv_dat$TRIALNO),]
britroc_surv_dat<-merge(surv_dat,pat_info,by.x=1,by.y=1)
rownames(britroc_surv_dat)<-britroc_surv_dat$IM.JBLAB_ID
britroc_surv_dat<-britroc_surv_dat[colnames(sig_pat_mat_britroc_all)[colnames(sig_pat_mat_britroc_all)%in%rownames(britroc_surv_dat)],]
britroc_surv_dat<-britroc_surv_dat[!(is.na(britroc_surv_dat$INT_START)|is.na(britroc_surv_dat$INT_END)),]
britroc_surv_dat<-britroc_surv_dat[!(britroc_surv_dat$INT_START>britroc_surv_dat$INT_END),]
britroc_train<-rownames(britroc_surv_dat)
britroc_test<-rownames(britroc_surv_dat)[!rownames(britroc_surv_dat)%in%britroc_train]
