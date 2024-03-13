#!/usr/bin/Rscript

# cmd <- paste("/misc/ngseq10/packages/Aligner/CellRanger/7.2.0/bin/cellranger mkref 
#              --genome Homo_sapiens_SMRTLinkgenome 
#              --fasta=/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/sl_unzip_datasets/5c640510-4bee-4ae8-b642-48ed209b79bb/call-unzip_datasets/execution/human_GRCh38_no_alt_analysis_set.fasta
#              --genes=/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/sl_unzip_datasets/5c640510-4bee-4ae8-b642-48ed209b79bb/call-unzip_datasets/execution/gencode.v39.annotation.sorted.gtf 
#              --memgb 20 
#              --output-dir /srv/GT/reference/Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_39-2021-12-09/Genes/genes_10XGEX_SC_SMRTLink_Index")
# ezSystem(cmd)

#Run CellRanger
library(stringr)
library(ezRun)

args = commandArgs(trailingOnly=TRUE)
print(args)
order = as.integer(args[1])
rawdirs = paste0("/srv/gstore/projects/p28443/", dir(path = "/srv/gstore/projects/p28443/", pattern = "o30669_NovaSeq_"))
sampledirs = list.files(path = rawdirs, pattern = "tar", full.names = TRUE)
samples = str_remove(basename(sampledirs), ".tar")
outputdir = "/srv/GT/analysis/zajacn/p28443/CellRanger_NoIntrons_SMRTLinkGTF"
refDir = "/srv/GT/reference/Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_39-2021-12-09/Genes/genes_10XGEX_SC_SMRTLink_Index"

setwd(outputdir)
for (i in sampledirs[order]){
  res <- untar(i, exdir = outputdir, tar=system("which tar", intern=TRUE))
  sample = str_remove(basename(i), ".tar")
  cmd <- paste(
    "/misc/ngseq10/packages/Aligner/CellRanger/7.2.0/bin/cellranger count", paste0("--id=", paste0(sample,"-cellRanger")),
    paste0("--transcriptome=", refDir),
    paste0("--fastqs=", paste0(outputdir,"/",sample)),
    paste0("--sample=", sample),
    paste0("--localmem=60"),
    paste0("--localcores=8"),
    paste0("--chemistry=auto"),
    paste0("--include-introns=false")
  )
  ezSystem(cmd)
} 

