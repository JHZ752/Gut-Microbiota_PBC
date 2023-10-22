
library(TwoSampleMR) 
library(dplyr)
library(utils)
library(data.table)
library(MRInstruments)
library(MRPRESSO)
library(RadialMR)
library(readxl) 
library(gridExtra) 
library(scales) 
library(tidyverse) 
library(tidyr)
library(withr)
index <- read.table("1.txt",as.is = T)#ieu ID
ID<-read.csv("microbiota.csv")
index
index <- as.vector(t(index))
all<- data.frame()
for (ebi in index) {exposure_dat <-extract_instruments(ebi,p1=1e-05,
                                                       clump = TRUE, r2 = 0.001,
                                                       kb = 10000,access_token = NULL)
outcome_dat <-extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST003129",
                                   access_token = NULL)
dat <- harmonise_data(exposure_dat,outcome_dat)
dat <- subset(dat,mr_keep)
res <- mr(dat,method_list = c("mr_ivw"))
res <-generate_odds_ratios(res)
res <-merge(ID,res,by="id.exposure")
#IVW
res$IVW.P<-res$pval
res$IVW.or<-res$or
res$IVW.or_lci95<-res$or_lci95
res$IVW.or_uci95<-res$or_uci95
res <- res%>% dplyr::select(id.exposure, Classification,Microbiota,ID,nsnp, IVW.P,IVW.or,IVW.or_lci95,IVW.or_uci95)
#egger
res_egger <- mr(dat,method_list = c("mr_egger_regression"))
res_egger <-generate_odds_ratios(res_egger)
res$egger.P<-res_egger$pval
res$egger.or<-res_egger$or
res$egger.or_lci95<-res_egger$or_lci95
res$egger.or_uci95<-res_egger$or_uci95
#WM
res_WM <- mr(dat,method_list = c("mr_weighted_median"))
res_WM <-generate_odds_ratios(res_WM)
res$WM.P<-res_WM$pval
res$WM.or<-res_WM$or
res$WM.or_lci95<-res_WM$or_lci95
res$WM.or_uci95<-res_WM$or_uci95
#MLE
res_MLE <- mr(dat,method_list = c("mr_two_sample_ml"))
res_MLE <-generate_odds_ratios(res_MLE)
res$MLE.P<-res_MLE$pval
res$MLE.or<-res_MLE$or
res$MLE.or_lci95<-res_MLE$or_lci95
res$MLE.or_uci95<-res_MLE$or_uci95
#WMODE
res_WMODE <- mr(dat,method_list = c("mr_weighted_mode"))
res_WMODE <-generate_odds_ratios(res_WMODE)
res$WMODE.P<-res_WMODE$pval
res$WMODE.or<-res_WMODE$or
res$WMODE.or_lci95<-res_WMODE$or_lci95
res$WMODE.or_uci95<-res_WMODE$or_uci95
####Sensitivity analysis####
pleio <- mr_pleiotropy_test(dat)
res$egger_intercept<-pleio$egger_intercept
res$egger_intercept_pval<-pleio$pval
###MRPRESSO results to confirm
results<- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                    SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                    data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
res$mrpresso.p <- results$`MR-PRESSO results`$`Global Test`$Pvalue
het <- mr_heterogeneity(dat)
het_egger<-het[1,]
het_IVW<-het[2,]
res$Q_pval_ivw<-het_IVW$Q_pval
res$Q_pval_mregger<-het_egger$Q_pval
table <- res%>% dplyr::select(id.exposure, Classification,Microbiota,ID,nsnp,
                              IVW.P,IVW.or,IVW.or_lci95,IVW.or_uci95,
                              egger.P,egger.or,egger.or_lci95,egger.or_uci95,
                              WM.P,WM.or,WM.or_lci95,WM.or_uci95,
                              MLE.P,MLE.or,MLE.or_lci95,MLE.or_uci95,
                              WMODE.P,WMODE.or,WMODE.or_lci95,WMODE.or_uci95,
                              egger_intercept,egger_intercept_pval,mrpresso.p,
                              Q_pval_ivw,Q_pval_mregger
                              )
if (nrow(res)!=0) {all <- rbind(all,table)
if(res$IVW.P<0.05){
  write.csv(table,file = paste0("",ebi,".csv"), row.names = FALSE)
}
}
}
write.csv(all,"all.csv",row.names=F)

