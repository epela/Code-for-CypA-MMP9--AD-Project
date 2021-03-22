rm(list = ls())

# Need devtools to install TwoSampleMR package 
# source("https://bioconductor.org/biocLite.R")
# to update MR-Base, run library(devtools), install_github and call the updated TWoSampleMR package
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR") #to update the package
devtools::install_github("MRCIEU/MRInstruments")
install.packages("psych")
install.packages("TwoSampleMR")
install.packages("MRInstruments")
install.packages("xlsx")
install.packages("ieugwasr")
install.packages("KNITR")
install.packages("xlsx")

library(ieugwasr)
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(xlsx)
library(psych)
library(TwoSampleMR)
library(ieugwasr)
library(knitr)

sessionInfo()

ao<-available_outcomes() 
#dim(ao)

########################################################################
########CYPA & MMP9 eQTL-->AD 2SMR ANALYSIS USING JANSEN AD GWAS########
########################################################################

setwd("C:/Users/epela/OneDrive - University of Bristol/MIGRATED O DRIVE/MRC Skills Development Fellowship/Projects/CYPA_MMP9/Data/eQTL analysis/2SMR")
#extracting eQTL instruments
CYPA_eQTLs <- extract_instruments("eqtl-a-ENSG00000196262")
write.csv(CYPA_eQTLs, "CypA_eQTLs.csv", row.names=T, quote = FALSE)
MMP9_eQTLs <- extract_instruments("eqtl-a-ENSG00000100985")
write.csv(MMP9_eQTLs, "MMP9_eQTLs.csv", row.names=T, quote = FALSE)

#reading in APOE SNP (association with AD but that's irrelevant as just need association with proteins)
apoe<- read_exposure_data("apoe_match.csv", sep = ",", snp_col="SNP", beta_col = "BETA", se_col="SE", effect_allele_col="A1", other_allele_col = "A2", eaf_col = "EAF")

#making CSV files for PRS generation in UKB
write.csv(CYPA_eQTLs, "CYPA_eQTLs_UKB.csv", row.names=T, quote = FALSE)
write.csv(MMP9_eQTLs, "MMP9_eQTLs_UKB.csv", row.names=T, quote = FALSE)


#reading in outcome data 
apoe_cypa_eQTL<- extract_outcome_data(apoe$SNP, c("eqtl-a-ENSG00000196262"), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
apoe_mmp9_eQTL<- extract_outcome_data(apoe$SNP, c("eqtl-a-ENSG00000100985"), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
CYPA_eQTLs_AD <- read_outcome_data("cypa_outcome_file.csv", sep = ",", snp_col="SNP", beta_col = "BETA", se_col="SE", effect_allele_col="A1", other_allele_col = "A2", eaf_col = "EAF")
MMP9_eQTLs_AD <- read_outcome_data("mmp9_outcome_file.csv", sep = ",", snp_col="SNP", beta_col = "BETA", se_col="SE", effect_allele_col="A1", other_allele_col = "A2", eaf_col = "EAF")


###################################
####EFFECT OF APOE ON CYPA & MMP9##
###################################
##READ OFF ASSOCIATION BETWEEN APOE & OUTCOME (CYPA/MMP9)
dat <- harmonise_data(apoe, apoe_cypa_eQTL)
dat2 <- harmonise_data(apoe, apoe_mmp9_eQTL)


###############################
###CYPA eQTL --> AD analysis###
###############################
dat <- harmonise_data(CYPA_eQTLs, CYPA_eQTLs_AD)

mr_results <- mr(dat)
mr_results
mr_report(dat, output_path = "CYPA_eQTLs_AD.TXT", author="Analyst", study = paste("CYPA_eQTLs","- Alzheimer's disease",sep=""))
results<-cbind.data.frame(mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)

#########################
###MMP9 eQTLs analysis###
#########################
dat2 <- harmonise_data(MMP9_eQTLs, MMP9_eQTLs_AD)

mr_results <- mr(dat2)
mr_results
mr_report(dat2, output_path = "MMP9_eQTLs_AD.txt", author="Analyst", study = paste("MMP9","- Alzheimer's disease",sep=""))
results<-cbind.data.frame(mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)


########################################################################
########CYPA & MMP9 pQTL-->AD 2SMR ANALYSIS USING JANSEN AD GWAS########
########################################################################
#NOTE: for APOE effect on CYPA and MMP9 pQTLs, look up APOE4 SNP in pQTLs SomaScan GWAS

setwd("C:/Users/epela/OneDrive - University of Bristol/MIGRATED O DRIVE/MRC Skills Development Fellowship/Projects/CYPA_MMP9/Data/pQTL analysis/2SMR")

#extracting pQTL instruments
CYPA_pQTLs <- read_exposure_data("GWS_SNPS_CYPA.csv", sep = ",", snp_col="snp", beta_col = "beta", se_col="se", effect_allele_col="effect_allele", other_allele_col = "other_allele", eaf_col = "effect_allele_freq")
CYPA_pQTLs<- clump_data(CYPA_exposure)

MMP9_pQTLs <- read_exposure_data("GWS_SNPS_MMP9.csv", sep = ",", snp_col="snp", beta_col = "beta", se_col="se", effect_allele_col="effect_allele", other_allele_col = "other_allele", eaf_col = "effect_allele_freq")
MMP9_pQTLs<- clump_data(MMP9_exposure)


#making CSV files for PRS generation in UKB
write.csv(CYPA_pQTLs, "CYPA_pQTLs_UKB.csv", row.names=T, quote = FALSE)
write.csv(MMP9_pQTLs, "MMP9_pQTLs_UKB.csv", row.names=T, quote = FALSE)


#reading in outcome data 
CYPA_pQTLs_AD <- read_outcome_data("cypa_pqtl_jansen_match.csv", sep = ",", snp_col="SNP", beta_col = "BETA", se_col="SE", effect_allele_col="A1", other_allele_col = "A2", eaf_col = "EAF")
MMP9_pQTLs_AD <- read_outcome_data("mmp9_pqtl_jansen_match.csv", sep = ",", snp_col="SNP", beta_col = "BETA", se_col="SE", effect_allele_col="A1", other_allele_col = "A2", eaf_col = "EAF")


###############################
###CYPA pQTL --> AD analysis###
###############################
dat <- harmonise_data(CYPA_pQTLs, CYPA_pQTLs_AD)

mr_results <- mr(dat)
mr_results
mr_report(dat, output_path = "CYPA_AD.txt", author="Analyst", study = paste("CYPA","- Alzheimer's disease",sep=""))
results<-cbind.data.frame(mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)

###############################
###MMP9 pQTL --> AD analysis###
###############################
dat2 <- harmonise_data(MMP9_pQTLs, MMP9_pQTLs_AD)

mr_results <- mr(dat2)
mr_results
mr_report(dat2, output_path = "MMP9_AD.txt", author="Analyst", study = paste("MMP9","- Alzheimer's disease",sep=""))
results<-cbind.data.frame(mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)
