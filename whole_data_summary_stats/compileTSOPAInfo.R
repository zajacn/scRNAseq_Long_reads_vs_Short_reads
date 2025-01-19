#!/usr/bin/Rscript
options(error=recover)
library(ezRun)
library(tidyverse)
library(parallel)
library(dplyr)
library(Rsamtools)
library(GenomicAlignments)
library(Seurat)

sample_data = read_tsv("../2024-05-22-sampleInfo.tsv")
sample_data = sample_data[5,]
workDir <- "/srv/GT/analysis/zajacn/p28443/"

compileTSOestimation_20240726 <- function(sampleInfo = sample_data, scratchDir="/srv/GT/analysis/zajacn/p28443/UMI_summary"){
  
  if (!file.exists(scratchDir)){
    stop()
  }
  sampleInfo$cellStatsFile <- paste0(scratchDir, "/", sampleInfo$Name, "-tso_pa.qsd")
  for (i in seq_along(sampleInfo$Name)){
    smi <- sampleInfo[i, ]
    message(smi$Name)
    if (!file.exists(smi$cellStatsFile)){
      illDir <- paste0("/srv/gstore/projects/", smi$CellRangerNoIntron)
      sbpill <- ScanBamParam(what =  c("qname", "rname"), tag = c("CB", "UB", "ts", "pa"))
      sbppb <- ScanBamParam(what =  c("qname", "rname"), tag = c("CB", "XM"), tagFilter=list(rc=1))
      
      pbBamList <- scanBam(smi$IsoSeqBam, param = sbppb)[[1]]
      pbumi = data.frame(pbBamList$tag)
      pbumi$CB = as.character(reverseComplement(DNAStringSet(pbumi$CB)))
      pbumi$XM = as.character(reverseComplement(DNAStringSet(pbumi$XM)))
      pbumi$CB = paste0(pbumi$CB, "-1")
      pbumi$tagID = paste(pbumi$CB, pbumi$XM)
      
      cts <- Read10X(file.path("/srv/gstore/projects", smi$CellRangerNoIntron, "filtered_feature_bc_matrix/"), gene.column = 1)
      illBamList <- scanBam(paste0(illDir, "/possorted_genome_bam.bam"), param = sbpill)[[1]]
      illumi = data.frame("qname" = illBamList$qname, illBamList$tag)
      illumi = illumi[illumi$CB %in% cts@Dimnames[[2]],]
      illumi = illumi[!duplicated(illumi),]
      illumi$tagID = paste(illumi$CB, illumi$UB)
      illumi$ts = if_else(illumi$ts == "NA", 0, 1)
      illumi$pa = if_else(illumi$pa == "NA", 0, 1)
      
      illumi$InPB = if_else(illumi$tagID %in% unique(pbumi$tagID), "InPB", "Unique")
      
      illumi = illumi %>% group_by(InPB, ts, pa) %>% dplyr::summarise(count = n_distinct(qname))
      qs::qsave(illumi, smi$cellStatsFile)
    }
    return(sampleInfo)
  }
}

sampleInfo <- compileTSOestimation_20240726(sampleInfo = sample_data, 
                                                 scratchDir="/srv/GT/analysis/zajacn/p28443/UMI_summary/")
