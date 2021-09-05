## ----setup, include=FALSE-----------------------------------------------------------------------------------------
source(".Rprofile")
#knitr::opts_chunk$set(echo = TRUE)
library(readr)
library(tidyverse)
library(viridis)
library(here)


## -----------------------------------------------------------------------------------------------------------------
eqtlgen <- read_tsv(here("data/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"))


## -----------------------------------------------------------------------------------------------------------------
#res_merge <- read_tsv("~/Documents/MPRA/Cancer_MPRA_Scaleup/AnalysisFiles12.5.19/DataAnalysis/Analysis_9.21.20/res_merged_9.21.20.tsv")
res_merge <- read_tsv(here("output/res_merge_hichip_nearbygenes_update_113020.tsv"))
# res_hmec <- read_tsv("~/Documents/MPRA/Cancer_MPRA_Scaleup/Sequencingv3/NextSeq_Pool1/Temp_HMEC_Analysis/mpranalyze_hmecNextSeq_oldpGF.tsv")
# res_merge <- right_join(res_hmec, res_merge)
res_merge <- dplyr::filter(res_merge, !is.na(Causal_SNP))


## -----------------------------------------------------------------------------------------------------------------
#sum(unique(eqtlgen$SNP) %in% res_merge$Causal_SNP)
#sum(unique(res_merge$Causal_SNP) %in% eqtlgen$SNP)



## -----------------------------------------------------------------------------------------------------------------
eg_snps <- distinct(dplyr::select(eqtlgen, SNP, SNPChr, SNPPos, AssessedAllele, OtherAllele))
# how many snps would be gained by using the position instead of rsid?
res_snps <- distinct(dplyr::select(res_merge, Causal_SNP, Chr, Start, Ref, Alt))

check_pos <- left_join(res_snps, eg_snps, by=c("Chr"="SNPChr", "Start"="SNPPos"))
sum(!is.na(check_pos$SNP))
# 3496
# 3828?
check_pos_id <- left_join(check_pos, dplyr::select(eg_snps, -AssessedAllele, -OtherAllele),
                          SNP, by=c("Causal_SNP"="SNP"))
sum(!is.na(check_pos_id$SNPPos))
# 3058
# 3083?
test <- check_pos_id[is.na(check_pos_id$SNPPos) & !is.na(check_pos_id$SNP),]
# using position saves some rsid inconsistencies, so add them as well


## -----------------------------------------------------------------------------------------------------------------
snpidlist <- unique(c(check_pos$Causal_SNP, check_pos$SNP))
# filter to only eatls in our library
eqtlgen <- dplyr::filter(eqtlgen, SNP %in% snpidlist)
# add a column to eqtlgen with the rsid the lib uses
eqtlgen <- eqtlgen %>% rowwise() %>%
  mutate(SNP_ID_lib = ifelse(SNP %in% check_pos$SNP, 
                             check_pos$Causal_SNP[check_pos$SNP == SNP & !is.na(check_pos$SNP)],
                             NA))
                             
            



## -----------------------------------------------------------------------------------------------------------------
check_alleles <- check_pos %>% rowwise() %>%
  mutate(eqtlgen_match = case_when(
    AssessedAllele == Ref & OtherAllele == Alt ~ "same",
    OtherAllele == Ref & AssessedAllele == Alt ~ "flip",
    is.na(AssessedAllele) ~ "no_eqtlgen",
    T ~ "eqtlgen_allele_mismatch"
  ))

# what information do I want on the eqtls?
# is there an eqtl? - eqtlgen_match, merge in on CausalSNP, Ref, Alt
# what are the egenes? - eqtlgen_egenes, collapse, semicolon separate
# what's the direction? - eqtlgen_dir, collapse, semicolon
# all same direction or diff? - eqtlgen_dir_sum, have same, diff




## -----------------------------------------------------------------------------------------------------------------
check_mult <- eqtlgen %>% group_by(SNP) %>% summarise(num_alleles = length(unique(AssessedAllele)))



## -----------------------------------------------------------------------------------------------------------------
check_alleles_onlymatch <- dplyr::filter(check_alleles, eqtlgen_match == "flip" | eqtlgen_match == "same")
eqtlgen_flip <- left_join(eqtlgen, dplyr::select(check_alleles_onlymatch, Causal_SNP, eqtlgen_match),
                          by=c("SNP_ID_lib"="Causal_SNP"))
collapse_eg_flip <- eqtlgen_flip %>% 
  mutate(Zscore_flip = ifelse(eqtlgen_match=="flip", -1*Zscore, Zscore),
         Ref = ifelse(eqtlgen_match == "flip", OtherAllele, AssessedAllele),
         Alt = ifelse(eqtlgen_match == "flip", AssessedAllele, OtherAllele)) %>%
  group_by(SNP_ID_lib, AssessedAllele, OtherAllele, Ref, Alt, eqtlgen_match) %>%
  summarise(eqtlgen_egenes = paste(GeneSymbol, collapse=";"),
            eqtlgen_dir_flipped =  paste(Zscore_flip, collapse=";"))



collapse_eg <- eqtlgen %>% group_by(SNP_ID_lib, AssessedAllele, OtherAllele) %>%
  summarise(eqtlgen_egenes = paste(GeneSymbol, collapse=";"),
            eqtlgen_dir_unflipped = paste(Zscore, collapse=";"))


## -----------------------------------------------------------------------------------------------------------------
res_merge_eg <- left_join(res_merge, 
                          dplyr::select(ungroup(collapse_eg_flip), SNP_ID_lib, Ref, Alt,
                                        eqtlgen_egenes, eqtlgen_dir_flipped),
                          by=c("Causal_SNP"="SNP_ID_lib", "Ref"="Ref", "Alt"="Alt"))

write_tsv(res_merge_eg, here("output/res_merge_witheQTLgen_11.1.2020.tsv"))
# res_merge_eg <- read_tsv(here("output/res_merge_witheQTLgen_11.1.2020.tsv"))

