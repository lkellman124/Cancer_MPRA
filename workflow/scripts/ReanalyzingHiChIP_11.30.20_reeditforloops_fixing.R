## ----setup, include=FALSE-----------------------------------------------------------------------------------------
source(".Rprofile")
#knitr::opts_chunk$set(echo = TRUE)
library(readr)
library(tidyverse)
library(here)


## -----------------------------------------------------------------------------------------------------------------
files = list.files(path = here("data/mergedloops/loops/"),
           full.names = T)
cell_types = str_remove(files, here("data/mergedloops/loops//"))
cell_types = str_remove(cell_types, "\\..*$")


## -----------------------------------------------------------------------------------------------------------------
# gene_table <- read_tsv("~/Documents/MPRA/Cancer_MPRA_Scaleup/eGeneFiles/HiChIP_112820/ucsc_gene_table.tsv")
gene_table <- read_tsv(here("data/ucsc_resfseq_refseqall"),
                       col_types = "ccccnncccccccccc")

# gene_conversion <- read_tsv("~/Documents/MPRA/Cancer_MPRA_Scaleup/eGeneFiles/HiChIP_112820/ucsc_id_conversion.tsv",
 #                           col_types = "cccccccccc")

# gene_table_merge <- left_join(gene_table, gene_conversion, by=c("#name"="#kgID"))
gene_table_merge <- gene_table
gene_table_merge <- dplyr::filter(gene_table_merge, !grepl("_", chrom))
gene_table_flip <- gene_table_merge %>% mutate(promoter = ifelse(strand=="+", txStart, txEnd))
  
gene_table_flip <- dplyr::select(gene_table_flip, name2, chrom, promoter) %>% distinct()

promoters <- gene_table_flip %>% group_by(name2, chrom) %>%
  summarise(prom_start = min(promoter) - 5000,
            prom_end = max(promoter) + 5000)
promoters <- promoters %>% mutate(chrom = str_remove(chrom, "chr"))
colnames(promoters)[1] <- "geneSymbol"


## -----------------------------------------------------------------------------------------------------------------
res_merge <- read_tsv(here("output/res_merge_plasmidrep_reanalyzed_112820.tsv"))
snp_positions <- dplyr::select(res_merge, Causal_SNP, Chr, Start) %>% distinct()
snp_positions$Start <- as.numeric(snp_positions$Start)


## -----------------------------------------------------------------------------------------------------------------
loops_all <- data.frame()
for (cell in cell_types){
  print(cell)
  loops_file <- read_csv(files[grep(cell, files)])
  loops_file <- loops_file %>% dplyr::select(-X1, -count) %>%
    separate(source, into=c("chr1", "s1", "e1"), remove=F) %>%
    separate(target, into=c("chr2", "s2", "e2"), remove=F)
  loops_file <- loops_file %>% mutate_at(c("s1", "e1", "s2", "e2"), as.numeric)
  loops_file <- loops_file %>% mutate(chr1 = str_remove(chr1, "chr"),
                                      chr2 = str_remove(chr2, "chr"))
  
  # Make a distinct anchor list
  anchors_1 <- distinct(dplyr::select(loops_file, source, chr1, s1, e1))
  colnames(anchors_1) <- c("anchor", "chrom", "start", "end")
  anchors_2 <- distinct(dplyr::select(loops_file, target, chr2, s2, e2))
  colnames(anchors_2) <- c("anchor", "chrom", "start", "end")
  anchors <- distinct(rbind(anchors_1, anchors_2))
  
  # Annotate anchor list with snps
  check_snp <- function(snp, snp_chrom, snp_pos, anchors) {
    anchor_snps <- dplyr::filter(anchors, (chrom == as.character(snp_chrom)) &
                                   (as.numeric(start) < as.numeric(snp_pos) & as.numeric(snp_pos) < as.numeric(end)))
    anchor_snps <- anchor_snps
    anchor_snps$snp <- unname(snp)
    anchor_snps$snp_pos <- unname(snp_pos)
    return(as.data.frame(anchor_snps))}
  
  
  check_snp <- function(snp, snp_chrom, snp_pos, anchors) {
    anchor_snps <- dplyr::filter(anchors, (chrom == as.character(snp_chrom)) &
                                   (as.numeric(start) < as.numeric(snp_pos) & as.numeric(snp_pos) < as.numeric(end)))
    anchor_snps <- anchor_snps %>% mutate(snp = snp, snp_pos = snp_pos)
    return(as.data.frame(anchor_snps))}
  
  # test <- check_snp("rs59351655", "8", 128351048, anchors)
  
  # check if apply is the problem
  snps_checked = data.frame()
  for (i in 1:dim(snp_positions)[1]){
    x = snp_positions[i,]
    snp_checked = check_snp(x["Causal_SNP"], x["Chr"], x["Start"], anchors)
    if (dim(snp_checked)[1] > 0){
      if(dim(snps_checked)[1] > 0){
        snps_checked = base::rbind(as.matrix(snps_checked), as.matrix(snp_checked))
      }else{
        snps_checked = snp_checked
      }
    }
  }
  anchor_snp_filt_df <- as.data.frame(snps_checked)
  
  anchor_snp_filt_df <- anchor_snp_filt_df %>% group_by(anchor, chrom, start, end) %>%
    summarise(snp = paste(unique(snp), collapse=";"))
  
  # Filter loop list down to those that touch snps
  loops_filt <- dplyr::filter(loops_file, source %in% anchor_snp_filt_df$anchor |
                                target %in% anchor_snp_filt_df$anchor)
  anchors_filt <- anchors[anchors$anchor %in% loops_filt$source |
                            anchors$anchor %in% loops_filt$target,]
  
  # change gene check to for loop too
  
  # Annotate anchor list with promoters
  check_prom <- function(gene, promoter_chrom, promoter_start, promoter_end, anchors) {
    anchor_genes <- dplyr::filter(anchors, (as.character(chrom) == as.character(promoter_chrom)) &
                                    ((as.numeric(start) < as.numeric(promoter_start) & 
                                        as.numeric(promoter_start) < as.numeric(end)) |
                                       (as.numeric(start) < as.numeric(promoter_end) & 
                                          as.numeric(promoter_end) < as.numeric(end)) |
                                       (as.numeric(promoter_start) < as.numeric(start) &
                                          as.numeric(end) < as.numeric(promoter_end))))
    if(dim(anchor_genes)[1] > 0){
      anchor_genes$gene <- gene
      anchor_genes <- anchor_genes %>% 
        mutate_at(c("anchor", "chrom"), as.character)
    }
    return(as.data.frame(anchor_genes))}
  
  genes_checked = data.frame()
  for (i in 1:dim(promoters)[1]){
    x = promoters[i,]
    gene_checked = check_prom(x[["geneSymbol"]], x[["chrom"]], x[["prom_start"]], x[["prom_end"]],
                              anchors_filt)
    if (dim(gene_checked)[1] > 0){
      if(dim(genes_checked)[1] > 0){
        genes_checked = base::rbind(as.matrix(genes_checked), as.matrix(gene_checked))
      }else{
        genes_checked = gene_checked
      }
    }
  }
  
  
  anchor_prom_filt_df <- as.data.frame(genes_checked)
  anchor_prom_filt_df <- anchor_prom_filt_df %>% group_by(anchor, chrom, start, end) %>%
    summarise(gene = paste(unique(gene), collapse=";"))
  
  # Combine the annotations
  anchors_annot <- full_join(anchor_snp_filt_df, anchor_prom_filt_df)
  
  # Filter to loops that have snps and genes
  loops_filt <- dplyr::filter(loops_filt, source %in% anchors_annot$anchor &
                                target %in% anchors_annot$anchor)
  loops_filt <- distinct(dplyr::select(loops_filt, source, target))
  loops_annot <- left_join(loops_filt, dplyr::select(ungroup(anchors_annot), anchor, snp, gene),
                           by = c("source"="anchor"))
  loops_annot <- left_join(loops_annot, dplyr::select(ungroup(anchors_annot), anchor, snp, gene),
                           by = c("target"="anchor"))
  loops_annot <- dplyr::filter(loops_annot, !(is.na(gene.x) & is.na(gene.y)))
  loops_annot$cell_type <- cell
  if(dim(loops_all)[1] > 0){
    loops_all = rbind(loops_all, loops_annot)
  } else{
    loops_all = loops_annot
  }
}


## -----------------------------------------------------------------------------------------------------------------
loops_all <- loops_all %>% 
  mutate(cell_type = case_when(
    grepl("GDS", cell_type) ~ "kc",
    grepl("Air", cell_type) ~ "airway",
    grepl("Ast", cell_type) ~ "ast",
    grepl("Colon", cell_type) ~ "colon",
    grepl("Uter", cell_type) ~ "endo",
    grepl("Eso", cell_type) ~ "eso",
    grepl("GM", cell_type) ~ "gm",
    grepl("HMEC", cell_type) ~ "hmec",
    grepl("Mel", cell_type) ~ "mc",
    grepl("Ov", cell_type) ~ "ov",
    grepl("Panc", cell_type) ~ "panc",
    grepl("Pros", cell_type) ~ "pros",
    grepl("Renal", cell_type) ~ "renal",
    grepl("Thy", cell_type) ~ "thy"))



## -----------------------------------------------------------------------------------------------------------------
# write_tsv(loops_all, here("output/loops_frommerged_112920.tsv"))
# write_tsv(loops_all, here("output/loops_rerunforloop_frommerged_113020.tsv"))
# write_tsv(loops_all, here("output/loops_rerunforloop_frommerged_113020_check.tsv"))
# loops_all_old <- read_tsv(here("output/loops_rerunforloop_frommerged_113020.tsv"))
# loops_all_old <- read_tsv(here("output/loops_frommerged_112920.tsv"))


## -----------------------------------------------------------------------------------------------------------------
loops_unnest <- loops_all %>% mutate(snp.x = str_split(snp.x, ";")) %>%
    unnest(cols=c(snp.x))
loops_unnest <- loops_unnest %>% mutate(snp.y = str_split(snp.y, ";")) %>%
    unnest(cols=c(snp.y))
loops_x = dplyr::select(loops_unnest, snp.x, gene.y, cell_type) %>%
  dplyr::filter(!is.na(snp.x)) %>% distinct()
loops_y = dplyr::select(loops_unnest, snp.y, gene.x, cell_type) %>%
  dplyr::filter(!is.na(snp.y)) %>% distinct()
colnames(loops_x) <- c("snp", "gene", "cell_type")
colnames(loops_y) <- c("snp", "gene", "cell_type")
loops <- rbind(loops_x, loops_y)
loops <- distinct(loops)
loops_snpcell <- loops %>% group_by(snp, cell_type) %>%
  summarise(gene = paste(unique(gene)[!is.na(unique(gene))], collapse=";"))
loops_snpcell <- dplyr::filter(loops_snpcell, !gene == "")


## -----------------------------------------------------------------------------------------------------------------
write_tsv(loops_snpcell, here("output/hichiploops_looponly_snpcellgrouped_113020.tsv"))
# loops_snpcell <- read_tsv("~/Documents/MPRA/Cancer_MPRA_Scaleup/eGeneFiles/HiChIP_112820/mergedloops/hichiploops_looponly_snpcellgrouped_112920.tsv")


## -----------------------------------------------------------------------------------------------------------------
same_bin_x = dplyr::select(loops_unnest, snp.x, gene.x, cell_type) %>%
  dplyr::filter(!is.na(snp.x) & !is.na(gene.x)) %>% distinct()
same_bin_y = dplyr::select(loops_unnest, snp.y, gene.y, cell_type) %>%
  dplyr::filter(!is.na(snp.y) & !is.na(gene.y)) %>% distinct()
colnames(same_bin_x) <- c("snp", "gene", "cell_type")
colnames(same_bin_y) <- c("snp", "gene", "cell_type")
same_bin <- rbind(same_bin_x, same_bin_y)
same_bin <- distinct(same_bin)
same_bin_snpcell <- same_bin %>% group_by(snp, cell_type) %>%
  summarise(gene = paste(unique(gene)[!is.na(unique(gene))], collapse=";"))
same_bin_snpcell <- dplyr::filter(same_bin_snpcell, !gene == "")



## -----------------------------------------------------------------------------------------------------------------
write_tsv(same_bin_snpcell, here("output/hichip_samebin_snpcellgrouped_112920.tsv"))
# same_bin_snpcell <- read_tsv("~/Documents/MPRA/Cancer_MPRA_Scaleup/eGeneFiles/HiChIP_112820/mergedloops/hichip_samebin_snpcellgrouped_112920.tsv")


## -----------------------------------------------------------------------------------------------------------------
snp_disease <- res_merge %>% dplyr::filter(!(Disease_expanded %in% c("bladder", "cervical"))) %>%
  dplyr::select(Causal_SNP, Disease_expanded) %>%
  mutate(Disease_expanded = str_split(Disease_expanded, ",")) %>%
  unnest(cols = Disease_expanded) %>% 
  mutate(Disease_expanded = case_when(
        grepl("kc", Disease_expanded) ~ "kc",
    grepl("airway", Disease_expanded) ~ "airway",
    grepl("ast", Disease_expanded) ~ "ast",
    grepl("colon", Disease_expanded) ~ "colon",
    grepl("endo", Disease_expanded) ~ "endo",
    grepl("eso", Disease_expanded) ~ "eso",
    grepl("lymphoma", Disease_expanded) ~ "gm",
    grepl("brca", Disease_expanded) ~ "hmec",
    grepl("mel", Disease_expanded) ~ "mc",
    grepl("ovca", Disease_expanded) ~ "ov",
    grepl("panc", Disease_expanded) ~ "panc",
    grepl("pros", Disease_expanded) ~ "pros",
    grepl("renal", Disease_expanded) ~ "renal",
    grepl("thy", Disease_expanded) ~ "thy"))

snp_disease <- snp_disease %>% group_by(Causal_SNP) %>%
  summarise(disease = paste(sort(unique(c(Disease_expanded))[!is.na(unique(c(Disease_expanded)))]), collapse=","))


## -----------------------------------------------------------------------------------------------------------------
loops_snpcell$loop <- T
same_bin_snpcell$loop <- F
# hichip <- rbind(loops_snpcell, same_bin_snpcell)
hichip <- loops_snpcell


## -----------------------------------------------------------------------------------------------------------------
hichip <- left_join(hichip, dplyr::select(snp_disease, Causal_SNP, disease), by=c("snp"="Causal_SNP"))
hichip <- dplyr::filter(hichip, !(disease==""))
hichip <- hichip %>% rowwise() %>% mutate(same_cell_type = grepl(cell_type, disease))


## -----------------------------------------------------------------------------------------------------------------
hichip_correct_celltype <- hichip %>% dplyr::filter(same_cell_type == T)
hichip_correct_celltype <- hichip_correct_celltype %>% group_by(snp) %>%
  summarise(hichip_correct_cell = paste(unique(gene), collapse=","),
            hichip_correct_cell_types = paste(unique(cell_type), collapse=","))
hichip_all_cell <- hichip %>% group_by(snp) %>%
  summarise(hichip_all = paste(gene, collapse=","),
            hichip_all_cell = paste(cell_type, collapse=","))


## -----------------------------------------------------------------------------------------------------------------
# write_tsv(hichip_correct_celltype, "~/Documents/MPRA/Cancer_MPRA_Scaleup/eGeneFiles/HiChIP_112820/mergedloops/hichiploops_looponly_correctcelltype_113020.tsv")
# hichip_correct_celltype_old <- read_tsv("~/Documents/MPRA/Cancer_MPRA_Scaleup/eGeneFiles/HiChIP_112820/mergedloops/hichiploops_looponly_correctcelltype_112920.tsv")
# write_tsv(hichip_all_cell, "~/Documents/MPRA/Cancer_MPRA_Scaleup/eGeneFiles/HiChIP_112820/mergedloops/hichiploops_looponly_allcells_113020.tsv")


## -----------------------------------------------------------------------------------------------------------------
res_merge_hichip <- left_join(res_merge, hichip_correct_celltype,
                              by=c("Causal_SNP"="snp"))
res_merge_hichip <- dplyr::select(res_merge_hichip, Causal_SNP,  
                                  Disease_expanded,hichip_correct_cell, hichip_correct_cell_types)
# test <- dplyr::filter(res_merge_hichip, !is.na(HiChIP_egene) | !is.na(hichip_egenes) |
#                         !is.na(hichip_correct_cell))


## -----------------------------------------------------------------------------------------------------------------
gene_loc <- dplyr::select(gene_table_merge, chrom, txStart, txEnd, name2)
gene_loc_gp <- gene_loc %>% group_by(name2) %>%
  summarise(chrom=unique(chrom),
            txStart = min(txStart),
            txEnd = max(txEnd))
# gene_loc <- gene_loc_gp %>% mutate(start = txStart - 10000,
#                                end = txEnd + 10000) 
gene_loc <- gene_loc_gp %>% mutate(start = txStart - 5000,
                                end = txEnd + 5000) 
gene_loc <- gene_loc %>% mutate(start = ifelse(start < 0, 0, start)) 
gene_loc <- gene_loc %>% mutate(chrom = str_remove(chrom, "chr"))
find_nearest <- function(snp, snp_chrom, snp_pos, gene_loc) {
  gene_snps <- dplyr::filter(gene_loc, (chrom == snp_chrom) &
                                (as.numeric(start) < as.numeric(snp_pos) & as.numeric(snp_pos) < as.numeric(end)))
  gene_snps$snp <- unname(snp)
  gene_snps$snp_pos <- unname(snp_pos)
  return(as.data.frame(gene_snps))}


nearest_genes_5k = data.frame()
for (i in 1:dim(snp_positions)[1]){
  x = snp_positions[i,]
  nearest_genes_snp = find_nearest(x[["Causal_SNP"]], x[["Chr"]], x[["Start"]], gene_loc)
  if (dim(nearest_genes_snp)[1] > 0){
      if(dim(nearest_genes_5k)[1] > 0){
    nearest_genes_5k = base::rbind(as.matrix(nearest_genes_5k), as.matrix(nearest_genes_snp))
  }else{
    nearest_genes_5k = nearest_genes_snp
  }
  }
}
nearest_genes <- as.data.frame(nearest_genes_5k)
nearest_genes_5k <- as.data.frame(nearest_genes_5k)

# using nearest_genes with 10kb filter
nearest_grouped <- nearest_genes %>% group_by(snp) %>%
  summarise(nearby_genes = paste(sort(unique(name2)), collapse=","))
  
  



## -----------------------------------------------------------------------------------------------------------------
colnames(hichip_correct_celltype) <- c("Causal_SNP", "hichip_egene", "hichip_cell_type")
res_merge_hichip <- left_join(res_merge,
                              hichip_correct_celltype)
colnames(nearest_grouped) <- c("Causal_SNP", "genes_within_10kb")
res_merge_gene <- left_join(res_merge_hichip, nearest_grouped)
check <-   dplyr::select(res_merge_gene, Causal_SNP, hichip_egene, genes_within_10kb)                                       
                                            


## -----------------------------------------------------------------------------------------------------------------
write_tsv(res_merge_gene, here("output/res_merge_hichip_nearbygenes_update_113020.tsv"))

