



# renv::load("~/Documents/MPRA/Cancer_MPRA_Scaleup/LK_CodeandDependencies/")

Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
# render the rmarkdowns

rmarkdown::render("Updating_merged_rsIDs.Rmd")
rmarkdown::render("HaploReg_FindingSNPs.Rmd")
rmarkdown::render("Finding_GWAS_RiskAllele_1.14.20.Rmd")
rmarkdown::render("GTex_VariantIDs.Rmd")
rmarkdown::render("GettingGTExGenes.Rmd")
rmarkdown::render("mpranalyze_plasmidrep_filter_11.27.20.Rmd")
rmarkdown::render("mpranalyze_merging_lib_background_info.Rmd")
rmarkdown::render("ReanalyzingHiChIP_11.30.20_reeditforloops_fixing.Rmd")
rmarkdown::render("IncorporatingeQTLGen.Rmd")
rmarkdown::render("alternatesharedegenevisualizations.Rmd")
rmarkdown::render("MarkingHitsandConvertingDiseases.Rmd")