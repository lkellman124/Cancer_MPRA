## ----setup, include=FALSE-----------------------------------------------------------------------------------------
source(".Rprofile")
#knitr::opts_chunk$set(echo = TRUE)
library(here)
library(tidyverse)


## -----------------------------------------------------------------------------------------------------------------
#lib <- read_tsv(here("output/lib_table_merged_rsids_12.10.19.tsv"))
 lib <- read_tsv(here("output/cancer_lib_haploreg_diseasereannotated_12.12.19.tsv"))
# lib_bg <- read.table(here("data/all-causal-SNP-16cancers.txt"), 
#                        stringsAsFactors = F, sep="\t", header=T )
gwas_cat <- read_tsv(here("data/gwas_catalog_v1.0-associations_e96_r2019-10-14.tsv"))


## -----------------------------------------------------------------------------------------------------------------
sum(grepl("\\?", gwas_cat$`STRONGEST SNP-RISK ALLELE`))
causal_snps <- unique(unlist(str_split(lib$Causal_SNP, ",")))
lead_snps <- unique(unlist(str_split(lib$leadSNP, ",")))

lib_snps <- unique(c(causal_snps, lead_snps))
gwas_lib_cat <- dplyr::filter(gwas_cat, SNPS %in% lib_snps)
                              
# checking to see if everything worked out                             
sum(lead_snps %in% gwas_lib_cat$SNPS)
sum(causal_snps %in% gwas_lib_cat$SNPS)

# 22 lead snps left out. Largely from a prostate interaction study that formatted the snps as snp x snp
# a lot of causal snps have no lead snp
# check <- lead_snps[!(lead_snps %in% gwas_lib_cat$SNPS)]
# checklib <- lib[lib$leadSNP %in% check,]
# checklib <- lib_bg[lib_bg$LeadSNP %in% check,]
# checklib <- checklib[!is.na(checklib$LeadSNP),]
gwas_filt <- gwas_cat
# checkgwas <- gwas_cat[unique(unlist(sapply(unique(checklib$LeadSNP), function(x) grep(x, gwas_cat$SNPS)))),]
# sapply(unique(checklib$LeadSNP), function(x) sum(grepl(x, checkgwas$SNPS)))

# all of the missing lead snps are also missing the useful allele information
# will just label them unknown

# what other disease are our snps linked to?
# check_disease <- dplyr::filter(gwas_lib_cat, !(grepl("cancer|carcinoma", `DISEASE/TRAIT`) | grepl("cancer", STUDY)))
# sort by disease so we don't get the allele linked to height or something
gwas_lib_filt <- gwas_lib_cat %>% dplyr::filter(grepl("ancer|arcinoma|ymphoma|lioma|eukemia|elanoma|lioblastoma|adenoma", `DISEASE/TRAIT`) | grepl("ancer", STUDY))


## -----------------------------------------------------------------------------------------------------------------
gwas_lib_filt <- gwas_lib_filt %>% mutate(gwas_allele = str_remove(`STRONGEST SNP-RISK ALLELE`,
                                                                 "^.*-"))



## -----------------------------------------------------------------------------------------------------------------
gwas_lib_filt <- gwas_lib_filt %>% 
  mutate(risk_allele = case_when(
    gwas_allele == "?" ~ "unknown",
    is.na(`OR or BETA`) ~ "unknown",
    grepl("increas", `95% CI (TEXT)`) ~ gwas_allele,
    grepl("decreas", `95% CI (TEXT)`) ~ paste("flip", gwas_allele, sep="_"),
    `OR or BETA` <  1 ~ "unknown",
    `OR or BETA` > 1 ~ gwas_allele,
    T ~ "unknown"            
  ))

duplicates <- gwas_lib_filt[gwas_lib_filt$SNPS %in% gwas_lib_filt$SNPS[duplicated(gwas_lib_filt$SNPS)],]
gwas_lib_coll <- gwas_lib_filt %>% group_by(SNPS) %>% 
  summarise(gwas_risk_allele = paste(risk_allele, collapse=","))
lead_gwas <- dplyr::select(gwas_lib_filt, SNPS, risk_allele)
           


## -----------------------------------------------------------------------------------------------------------------
# write_tsv(gwas_lib_coll, here("output/gwas_risk_allele_key_leadsnps_veryunreliable_1.14.20.tsv"))


## -----------------------------------------------------------------------------------------------------------------
haplo <- read_tsv(here("data/cmpra_haploreg_data_withquery.tsv"))
rc_check <- function(x) {
  return(case_when(
    x == "A" ~ "T",
    x == "T" ~ "A",
    x == "C" ~ "G",
    x == "G" ~ "C"
  ))
}

haplo_gwas <- right_join(haplo, lead_gwas, by=c("rsID"="SNPS"))
haplo_gwas <- haplo_gwas %>% mutate(risk_allele_ref_alt = case_when(
  risk_allele == "unknown" ~ "unknown",
  risk_allele == ref ~ "ref",
  risk_allele == alt ~ "alt",
  !grepl(",", alt) & rc_check(risk_allele) == ref ~ "ref",
  !grepl(",", alt) & rc_check(risk_allele) == alt ~ "alt",
  T ~ "unknown"))

# write_tsv(haplo_gwas, here("output/lead_snp_risk_allele_haplo_gwas_comparisons_1.14.20.tsv"))
# haplo_gwas <- read_tsv("lead_snp_risk_allele_haplo_gwas_comparisons_1.14.20.tsv")


## -----------------------------------------------------------------------------------------------------------------
# lib <- read_tsv("~/Documents/MPRA/Cancer_MPRA_Scaleup/AnalysisFiles12.5.19/FormattingBackground/lib_studies_idmerge_diseasesexpanded_coords_GTEx_hichip_1.13.20_updated.tsv")

haplo_gwas_coll <- haplo_gwas %>% dplyr::select(rsID, risk_allele_ref_alt) %>% distinct() %>%
  group_by(rsID) %>% summarise(gwas_risk_allele = paste("", risk_allele_ref_alt, sep="", collapse=","))
lib_merge <- left_join(lib, haplo_gwas_coll, by=c("leadSNP"="rsID"))
write_tsv(lib_merge, here("output/lib_studies_idmerge_diseasesexpanded_coords_gtex_hichip_gwasrisk_1.14.20.tsv"))

