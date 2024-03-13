#!/usr/bin/Rscript
options(error=recover)
library(Seurat)
library(HDF5Array)
library(ezRun)
library(RColorBrewer)
library(SingleCellExperiment)
library(pheatmap)
library(DropletUtils)
library(scuttle)
library(BiocParallel)
library(data.table)
library(ggplot2)
library(ezRun)
library(Seurat)
library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)


samples <- ezRead.table("/home/zajac/Genomics/p28443_singlecell_pacbiovsillumina/2024-01-17-sampleInfo.tsv")
workDir <- "/srv/GT/analysis/zajacn/p28443/"

for (sm in rownames(samples)[2:8]){
  
  sampleInfo <- samples[sm, ] %>% as.list()
  rmarkdown::render("/srv/GT/analysis/zajacn/p28443/2024-06-01-bams-cell-subset-multimapping.Rmd", 
                    output_file = paste0(dirname(sampleInfo$PigeonUnfiltered), "/2024-01-29-multimapping-reads-", sm, ".html"), 
                    output_dir = dirname(sampleInfo$PigeonUnfiltered), knit_root_dir=workDir)
}