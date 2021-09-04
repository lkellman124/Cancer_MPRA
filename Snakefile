AREAS = "Brain Colon EBV Esop Thyroid Skin Pancreas Breast Lung Uterus Kidney Ovary Prostate".split()

rule all:
  input: "PLACEHOLDER"

rule make_res_ids:
  input: 
    "data/ref_alt_snp_table.txt",
    "data/RsMergeArch.bcp"
  output:
    "output/lib_snps_merged_rsid_conversion.tsv",
    "output/lib_table_merged_rsids_12.10.19.tsv"
  script:
    "Updating_merged_rsIDs.Rmd"
    
rule gtex_variant_ids:
  input:
    "data/lib_studies_idmerge_diseasesexpanded_coords3738_gtexID_12.13.19.tsv",
    "data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt"
  output:
    "output/lib_studies_idmerge_diseasesexpanded_coords3738_gtexIDuncollapsed_1.9.20.tsv"
  script:
    "GTex_VariantIDs.Rmd"
    
rule gtex_genes:
  input:
    "output/lib_studies_idmerge_diseasesexpanded_coords3738_gtexIDuncollapsed_1.9.20.tsv",
    expand("data/LK_GTEx_Analysis_v8_eQTL/{area}.v8.signif_variant_gene_pairs.txt", area=AREAS)
  output:
    "output/significant_egene_pairs_for_snps_in_gtex_tissue_filtered_1.13.20.txt",
    "output/lib_studies_idmerge_diseasesexpanded_coords3738_GTExeGenes_tissuecolumn_1.13.20.tsv"
  script:
    "GettingGTExGenes.Rmd"

rule mpra_analyze_plasmidrep:
  input:
    "data/nextseq_counts_111920.tsv"
  output:
    "output/mpranalyze_plasmidrep_112720.tsv"
  script:
    "mpranalyze_plasmidrep_filter_11.27.20.Rmd"

rule mpra_analyze_merge_background_info:
  input:
    "output/mpranalyze_plasmidrep_112720.tsv",
    "data/library_risk_allele_adjusted_locus_matchup.tsv"
  output:
    "output/res_merge_plasmidrep_reanalyzed_112820.tsv"
  script:
    "mpranalyze_merging_lib_background_info.Rmd"