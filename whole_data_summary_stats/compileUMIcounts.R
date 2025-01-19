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
sample_data$IsoSeqBam2 = list.files("/scratch/zajacn/p28443/all_cells/", pattern = "scisoseq.mapped.bam$", recursive = TRUE, full.names = TRUE) #consisting of all cells
sample_data = sample_data[5,]
workDir <- "/srv/GT/analysis/zajacn/p28443/"

compileUMIperCellrecovery_20240709 <- function(sampleInfo = sample_data, scratchDir="/srv/GT/analysis/zajacn/p28443/UMI_summary"){
   
   if (!file.exists(scratchDir)){
     stop()
   }
   i <- 1
   sampleInfo$cellStatsFile <- paste0(scratchDir, "/", sampleInfo$Name, "-cellUMIStats.qsd")
   for (i in seq_along(sampleInfo$Name)){
     smi <- sampleInfo[i, ]
     message(smi$Name)
     if (!file.exists(smi$cellStatsFile)){
       illDir <- paste0("/srv/gstore/projects/", smi$CellRangerNoIntron)
       sbpill <- ScanBamParam(what =  c("qname", "rname"), tag = c("CB", "UB"))
       sbppb <- ScanBamParam(what =  c("qname", "rname"), tag = c("CB", "XM"), tagFilter=list(rc=1))
       sbppb2 <- ScanBamParam(what =  c("qname", "rname"), tag = c("CB", "XM"))
       
       pbBamList <- scanBam(smi$IsoSeqBam, param = sbppb)[[1]] #consisting of real cells only
       pbumi = data.frame(pbBamList$tag)
       pbumi$CB = as.character(reverseComplement(DNAStringSet(pbumi$CB)))
       pbumi$CB = paste0(pbumi$CB, "-1")
       pbumi$TotalPB = length(unique(paste(pbumi$CB, pbumi$XM)))
       
       # For the cells missing from PacBio and for comparison with the bam file after retaining nonreal cells
       pbBamList <- scanBam(smi$IsoSeqBam2, param = sbppb2)[[1]]
       pbumi2 = data.frame(pbBamList$tag)
       pbumi2$CB = as.character(reverseComplement(DNAStringSet(pbumi2$CB)))
       pbumi2$CB = paste0(pbumi2$CB, "-1")
       pbumi2$TotalPB = length(unique(paste(pbumi2$CB, pbumi2$XM)))
       
       cts <- Read10X(file.path("/srv/gstore/projects", smi$CellRangerNoIntron, "filtered_feature_bc_matrix/"), gene.column = 1)
       illBamList <- scanBam(paste0(illDir, "/possorted_genome_bam.bam"), param = sbpill)[[1]]
       illumi = data.frame(illBamList$tag)
       illumi = illumi[illumi$CB %in% cts@Dimnames[[2]],]
       illumi = illumi[!duplicated(illumi),]
       illumi$TotalIll = length(unique(paste(illumi$CB, illumi$UB)))
       
       illumi$Commonness <- if_else(illumi$CB %in% pbumi$CB, "Common", "Unique to ILL")
       pbumi$Commonness <- if_else(pbumi$CB %in% illumi$CB, "Common", "Unique to PB")
       pbumi2  = pbumi2[pbumi2$CB %in% illumi$CB,]
       
       illumi = illumi %>% group_by(CB, Commonness, TotalIll) %>% dplyr::summarise(UMI_count_Ill = n())
       pbumi = pbumi %>% group_by(CB, Commonness, TotalPB) %>% dplyr::summarise(UMI_count_Pb_real = n())
       pbumi2 = pbumi2 %>% group_by(CB, TotalPB) %>% dplyr::summarise(UMI_count_Pb_all = n())
       
       xt = merge(pbumi, illumi, by = c("CB", "Commonness"), all = TRUE) %>% column_to_rownames("CB") %>% replace(is.na(.), 0)
       xt = merge(xt, pbumi2, by = "CB", all = TRUE)
       
       qs::qsave(xt, smi$cellStatsFile)
     }
   return(sampleInfo)
   }
}

sampleInfo <- compileUMIperCellrecovery_20240709(sampleInfo = sample_data, 
                                           scratchDir="/srv/GT/analysis/zajacn/p28443/UMI_summary/")