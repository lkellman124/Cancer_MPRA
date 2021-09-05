## ----setup, include=FALSE-----------------------------------------------------------------------------------------
source(".Rprofile")
#knitr::opts_chunk$set(echo = TRUE)
library(readr)
library(magrittr)
library(dplyr)
library(here)


## -----------------------------------------------------------------------------------------------------------------
# lib <- read_tsv("~/Documents/MPRA/Cancer_MPRA_Scaleup/AnalysisFiles12.5.19/FormattingBackground/GenomicCoordinates/lib_studies_idmerge_diseasesexpanded_coords3738_12.13.19.tsv")
lib <- read_tsv(here("output/lib_studies_idmerge_diseasesexpanded_coords_gtex_hichip_gwasrisk_1.14.20.tsv"))

gtex_lookup <- read_tsv(here("data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt"))




## -----------------------------------------------------------------------------------------------------------------
all_snps = unique(c(lib$Causal_SNP, lib$Lead_SNP))

gtex_rsids <- gtex_lookup %>% filter(rs_id_dbSNP151_GRCh38p7 %in% all_snps)
gtex_leftover <- all_snps[!(all_snps %in% gtex_rsids$rs_id_dbSNP151_GRCh38p7)]

# 426 total leftover
# 341 causal snps


## -----------------------------------------------------------------------------------------------------------------
gtex_rsids_causal <- dplyr::filter(gtex_rsids, rs_id_dbSNP151_GRCh38p7 %in% lib$Causal_SNP) %>%
  dplyr::select(variant_id, rs_id_dbSNP151_GRCh38p7, num_alt_per_site, alt) 

# need to deal with the duplicates in a way that won't interfere later with merging in e-genes.
# Maybe allow multiple entries

lib <- left_join(lib, gtex_rsids_causal,
                 by=c("Causal_SNP" = "rs_id_dbSNP151_GRCh38p7"))
# expands it by a few
write_tsv(lib, here("output/lib_studies_idmerge_diseasesexpanded_coords3738_gtexIDuncollapsed_1.9.20.tsv")) 
# gtex_rsids_causal_sum <- gtex_rsids_causal %>% group_by(rs_id_dbSNP151_GRCh38p7) %>%
#   summarise(gtex_alt = paste(unique(alt), collapse=","), 
#             gtex_num_alt_per_site = unique(num_alt_per_site),
#             gtex_variant_id = paste(sort(unique(variant_id)), collapse=","))

# lib <- left_join(lib, gtex_rsids_causal_sum,
                 # by=c("Causal_SNP" = "rs_id_dbSNP151_GRCh38p7"))
   

# write_tsv(lib, "~/Documents/MPRA/Cancer_MPRA_Scaleup/AnalysisFiles12.5.19/FormattingBackground/GTEx_IDs_eQTLs/lib_studies_idmerge_diseasesexpanded_coords3738_gtexID_12.13.19.tsv")

