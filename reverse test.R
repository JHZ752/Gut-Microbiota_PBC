#根据结局找暴露
##############加载包
library(TwoSampleMR) 
library(dplyr)
library(utils)
library(TwoSampleMR)
library(data.table)
library(readxl) 
library(gridExtra) 
library(scales) 
library(tidyverse) 
ID<-read.csv("microbiota.csv")
#exposure_dat <-extract_instruments("ebi-a-GCST003129",p1=1e-05,
                                   #clump = TRUE, r2 = 0.001,
                                   #kb = 10000,access_token = NULL)
#exposure_dat<-snp_add_eaf(exposure_dat)
#write.csv(exposure_dat,"PBCEXP.csv")
exposure_dat<-read.csv("PBCEXP.csv",header = T,row.names = 1)

outcome_dat <-extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST90016946",access_token = NULL)#选择结局
#去除中等等位基因频率的SNP
dat <- harmonise_data(exposure_dat,outcome_dat)
dat <- subset(dat,mr_keep)
#排除与结局相关的SNP
dat<- subset(dat,pval.outcome>=5e-08)
res <- mr(dat,method_list = c("mr_ivw"))#mr方法选择
dat$EAF2 <- (1 - dat$eaf.exposure)
dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)
PVEfx <- function(BETA, MAF, SE, N){
  pve <- (2*(BETA^2)*MAF*(1 - MAF))/((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF)))
  return(pve) 
}
dat$R2 <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.exposure)
dat$FSTAT <- ((dat$samplesize.exposure - 1 - 1)/1)*(dat$R2/(1 - dat$R2))
dat<-merge(ID,dat,by="id.outcome")
table <- dat%>% dplyr::select(id.outcome,Classification,Microbiota,ID,SNP,effect_allele.exposure, other_allele.exposure,beta.exposure,se.exposure,pval.exposure,eaf.exposure,
                              beta.outcome,se.outcome,pval.outcome,R2, FSTAT)
write.csv(table,"ebi-a-GCST90016946.csv")

