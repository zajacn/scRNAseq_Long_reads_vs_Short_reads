#!/usr/bin/Rscript
options(error=recover)
library(ezRun)
library(tidyverse)
library(parallel)
library(dplyr)
library(Rsamtools)
library(GenomicAlignments)
library(Seurat)

sample_data = read_tsv("~/Genomics/p28443_singlecell_pacbiovsillumina/SIB-poster/2024-05-22-sampleInfo.tsv")
sample_data = sample_data[5,]
workDir <- "/srv/GT/analysis/zajacn/p28443/"

compileCellStats_20240805 <- function(sampleInfo = sample_data,
                                      scratchDir="/srv/GT/analysis/zajacn/p28443/UMI_summary/"){
  
  if (!file.exists(scratchDir)){
    stop()
  }
  i <- 1
  sampleInfo$cellStatsFile <- paste0(scratchDir, "/", sampleInfo$Name, "-cellMITOStats.qsd")
  for (i in seq_along(sampleInfo$Name)){
    smi <- sampleInfo[i, ]
    message(smi$Name)
    if (!file.exists(smi$cellStatsFile)){
      illDir <- paste0("/srv/gstore/projects/", smi$CellRangerNoIntron)
      sbp <- ScanBamParam(what =  c("qname", "rname", "qwidth"), tag = c("CB"))
      pbBamList <- scanBam(smi$IsoSeqBam, param = sbp)[[1]]
      isMito <- pbBamList$rname %in% "chrM"
      pbPercentMito <- tapply(isMito, pbBamList$tag$CB, mean) * 100
      names(pbPercentMito) = as.character(reverseComplement(DNAStringSet(names(pbPercentMito))))
      
      cts <- Read10X(file.path("/srv/gstore/projects", smi$CellRangerNoIntron, "filtered_feature_bc_matrix/"), gene.column = 1)
      illBamList <- scanBam(paste0(illDir, "/possorted_genome_bam.bam"), param = sbp)[[1]]
      use <- illBamList$tag$CB %in% cts@Dimnames[[2]]
      xt <- table(illBamList$rname[use], illBamList$tag$CB[use], useNA="always")
      illPercentMito <- xt["chrM", ]/colSums(xt) * 100
      names(illPercentMito) <- sub("-1", "", names(illPercentMito))
      cellStats <- merge(as.data.frame(pbPercentMito) %>% rownames_to_column("CellID"), as.data.frame(illPercentMito) %>% rownames_to_column("CellID"), by = "CellID", all = TRUE)
      qs::qsave(cellStats, smi$cellStatsFile)
    }
  }
  return(sampleInfo)
}

sampleInfo <- compileCellStats_20240805(sampleInfo = sample_data, 
                                                 scratchDir="/srv/GT/analysis/zajacn/p28443/UMI_summary/")
