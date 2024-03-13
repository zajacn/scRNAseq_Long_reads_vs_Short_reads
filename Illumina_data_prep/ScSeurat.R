#!/usr/bin/Rscript
library(HDF5Array)
library(AUCell)
library(GSEABase)
library(SingleR)
library(Seurat)
library(tidyverse)
library(scDblFinder)
library(BiocParallel)
library(scuttle)
library(DropletUtils)
library(enrichR)
library(decoupleR)
library(Azimuth)
library(ezRun)

###This script runs seurat on Illumina data using the SMRTLink v11.1 annotation

##Create annotation files
# setwd("/srv/GT/reference/Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_39-2021-12-09/Genes/SMRTLink11.1_annotation/")
# genomeFile="/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/sl_unzip_datasets/5c640510-4bee-4ae8-b642-48ed209b79bb/call-unzip_datasets/execution/human_GRCh38_no_alt_analysis_set.fasta"
# featureFile="/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/sl_unzip_datasets/5c640510-4bee-4ae8-b642-48ed209b79bb/call-unzip_datasets/execution/gencode.v39.annotation.sorted.gtf"
# 
# 
# makeFeatAnnoEnsemblmodified <- function(featureFile,
#          genomeFile,
#          biomartFile=NULL,
#          organism=NULL,
#          host=NULL,
#          mart='ENSEMBL_MART_ENSEMBL'){
#   require(rtracklayer)
#   require(data.table)
#   
#   featAnnoFile <-  paste0("genes_annotation_byTranscript.txt")
#   featAnnoGeneFile <- paste0("genes_annotation_byGene.txt")
#   
#   feature <- import(featureFile)
#   transcripts <- feature[feature$type=="transcript"]
#   if(length(transcripts) == 0L){
#     ## Incomplete gtf with only exons.
#     ## Try to reconstruct the transcripts.
#     exons <- feature[feature$type == "exon"]
#     exonsByTx <- GenomicRanges::split(exons, exons$transcript_id)
#     transcripts <- unlist(range(exonsByTx))
#     transcripts$transcript_id <- names(transcripts)
#     names(transcripts) <- NULL
#     transcripts$gene_id <- exons$gene_id[match(transcripts$transcript_id, 
#                                                exons$transcript_id)]
#     transcripts$gene_name <- exons$gene_name[match(transcripts$transcript_id, 
#                                                    exons$transcript_id)]
#     transcripts$gene_type <- exons$gene_type[match(transcripts$transcript_id, 
#                                                          exons$transcript_id)]
#   }
#   
#   transcripts <- transcripts[!duplicated(transcripts$transcript_id)]
#   ## This is to deal with the cases of duplicates transcripts from GENCODE annotation
#   ## Example: ENST00000399012 can be on chrX and chrY.
#   ## Ensembl only keeps the ones on chrX.
#   
#   ## Calculate gc and featWidth
#   gw <- getTranscriptGcAndWidth(genomeFn=genomeFile,
#                                 featureFn=featureFile)
#   gw$transcript_id = str_replace(gw$transcript_id, pattern = ".[0-9]+$",replacement = "")
#   featAnno <- tibble(transcript_id=str_replace(transcripts$transcript_id, pattern = ".[0-9]+$",replacement = ""),
#                      gene_id=str_replace(transcripts$gene_id, pattern = ".[0-9]+$",replacement = ""),
#                      gene_name=transcripts$gene_name,
#                      type=transcripts$gene_type,
#                      strand=as.character(strand(transcripts)),
#                      seqid=as.character(seqnames(transcripts)),
#                      start=start(transcripts),
#                      end=end(transcripts),
#                      biotypes=transcripts$gene_type)
#   featAnno <- left_join(featAnno, gw)
#   
#   ## The numeric columns should not have NAs
#   stopifnot(!featAnno %>% dplyr::select(start, end, gc, featWidth) %>% 
#               is.na() %>% any())
#   
#   ## Group the biotype into more general groups
#   stopifnot(all(featAnno %>% pull(biotypes) %in% listBiotypes("all")))
#   isProteinCoding <- featAnno %>% pull(biotypes) %in% listBiotypes("protein_coding")
#   isLNC <- featAnno %>% pull(biotypes) %in% listBiotypes("long_noncoding")
#   isSHNC <- featAnno %>% pull(biotypes) %in% listBiotypes("short_noncoding")
#   isrRNA <- featAnno %>% pull(biotypes) %in% listBiotypes("rRNA")
#   istRNA <- featAnno %>% pull(biotypes) %in% listBiotypes("tRNA")
#   isMtrRNA <- featAnno %>% pull(biotypes) %in% listBiotypes("Mt_rRNA")
#   isMttRNA <- featAnno %>% pull(biotypes) %in% listBiotypes("Mt_tRNA")
#   isPseudo <- featAnno %>% pull(biotypes) %in% listBiotypes("pseudogene")
#   featAnno$type[isPseudo] <- "pseudogene"
#   featAnno$type[isLNC] <- "long_noncoding"
#   featAnno$type[isSHNC] <- "short_noncoding"
#   featAnno$type[isProteinCoding] <- "protein_coding"
#   ### rRNA and tRNA have to be after noncoding
#   ### since they are subset of noncoding
#   featAnno$type[isrRNA] <- "rRNA"
#   featAnno$type[istRNA] <- "tRNA"
#   featAnno$type[isMtrRNA] <- "Mt_rRNA"
#   featAnno$type[isMttRNA] <- "Mt_tRNA"
#   
#   ## additional information from Ensembl or downloaded biomart file
#   attributes <- c("ensembl_transcript_id", "description", 
#                   "go_id", "namespace_1003")
#   names(attributes) <- c("Transcript stable ID", "Gene description",
#                          "GO term accession", "GO domain")
#   ## Older web-page biomart has different names
#   attributesOld <- set_names(attributes,
#                              c("Ensembl Transcript ID", "Description",
#                                "GO Term Accession", "GO domain"))
#   if(!is.null(biomartFile)){
#     message("Using local biomart file!")
#     # fread cannot handle compressed file
#     mapping <- as.data.table(read_tsv(biomartFile, guess_max=1e6)) 
#     if(all(names(attributes) %in% colnames(mapping))){
#       mapping <- mapping[ ,names(attributes), with=FALSE]
#       # To make it consistent with biomaRt
#       colnames(mapping) <- attributes[colnames(mapping)] 
#     }else if(all(names(attributesOld) %in% colnames(mapping))){
#       mapping <- mapping[ ,names(attributesOld), with=FALSE]
#       # To make it consistent with biomaRt
#       colnames(mapping) <- attributesOld[colnames(mapping)]
#     }else{
#       stop("Make sure ", paste(names(attributes), collapse="; "), 
#            "are downloaded from web biomart!")
#     }
#   }else if(!is.null(organism)){
#     message("Query via biomaRt package!")
#     require(biomaRt)
#     if(is.null(host)){
#       ensembl <- useMart(mart)
#     }else{
#       ensembl <- useMart(mart, host=host)
#     }
#     ensembl <- useDataset(organism, mart=ensembl)
#     mapping1 <-
#       getBM(attributes=setdiff(attributes, c("go_id", "namespace_1003")),
#             filters=c("ensembl_transcript_id"),
#             values=featAnno$transcript_id, mart=ensembl)
#     mapping1 <- as_tibble(mapping1)
#     mapping2 <-
#       getBM(attributes=c("ensembl_transcript_id", "go_id", "namespace_1003"),
#             filters=c("ensembl_transcript_id"),
#             values=featAnno$transcript_id, mart=ensembl)
#     mapping2 <- as_tibble(mapping2)
#     mapping <- inner_join(mapping1, mapping2)
#   }else{
#     message("Not using any additional annotation!")
#     mapping <- tibble(ensembl_transcript_id=featAnno$transcript_id,
#                       description="", go_id="", namespace_1003="")
#   }
#   mapping <- mapping %>%
#     mutate(ensembl_transcript_id=str_replace(ensembl_transcript_id, "\\.\\d+$", ""))
#   
#   if(!all(featAnno$transcript_id %in% mapping$ensembl_transcript_id)){
#     warning("Some transcript ids don't exist in biomart file!") #Normal for GENCODE
#   }
#   
#   ### description
#   txid2description <- mapping %>% dplyr::select(transcript_id=ensembl_transcript_id,
#                                                 description) %>%
#     filter(!duplicated(transcript_id))
#   featAnno <- left_join(featAnno, txid2description)
#   
#   ### GO
#   GOMapping <- c("biological_process"="GO BP",
#                  "molecular_function"="GO MF",
#                  "cellular_component"="GO CC")
#   go <- mapping %>% dplyr::select(transcript_id=ensembl_transcript_id, go_id, namespace_1003) %>%
#     filter(!(is.na(go_id) | go_id == ""),
#            namespace_1003 %in% c("biological_process", "molecular_function", "cellular_component")) %>%
#     group_by(transcript_id, namespace_1003) %>%
#     summarise(go_id=str_c(unique(go_id), collapse="; ")) %>% ungroup()
#   if(nrow(go)==0L){
#     ## If go is an empty data.table
#     go <- tibble(transcript_id=featAnno$transcript_id,
#                  biological_process="", molecular_function="", cellular_component="")
#   }else{
#     go <- pivot_wider(go, id_cols=transcript_id, names_from = namespace_1003,
#                       values_from=go_id, values_fill="")
#   }
#   go <- dplyr::rename(go, "GO BP"="biological_process", "GO MF"="molecular_function",
#                       "GO CC"="cellular_component")
#   featAnno <- left_join(featAnno, go)
#   featAnno <- featAnno %>% mutate("GO BP"=replace_na(`GO BP`, ""),
#                                   "GO MF"=replace_na(`GO MF`, ""),
#                                   "GO CC"=replace_na(`GO CC`, ""))
#   
#   ## output annotation file on transcript level
#   write_tsv(featAnno, file=featAnnoFile)
#   
#   ## make annotation at gene level
#   featAnnoGene <- aggregateFeatAnno(featAnno)
#   write_tsv(featAnnoGene, file=featAnnoGeneFile)
#   
#   invisible(list("transcript"=featAnno, "gene"=featAnnoGene))
# }
# 
# makeFeatAnnoEnsemblmodified(featureFile=featureFile, genomeFile=genomeFile, organism="hsapiens_gene_ensembl", host="nov2020.archive.ensembl.org")

#Define functions (copy from ScSeurat App in SUSHI)
addCellQcToSeurat <- function(scData, param=NULL, BPPARAM=NULL, ribosomalGenes=NULL){
  
  library(scater)
  
  # Cells filtering
  scData <- PercentageFeatureSet(scData, "(?i)^MT-", col.name = "percent_mito")
  scData <- PercentageFeatureSet(scData, "(?i)^RPS|^RPL", col.name = "percent_riboprot")
  if (!is.null(ribosomalGenes)){
    scData <- PercentageFeatureSet(scData, features=ribosomalGenes, col.name = "percent_ribosomal")
  }
  if(grepl("Spatial", param$appName)) {
    assay <- "Spatial"
    att_nCounts <- "nCount_Spatial"
    att_nGenes <- "nFeature_Spatial"
  } else {
    att_nCounts <- "nCount_RNA"
    att_nGenes <- "nFeature_RNA"
    assay <- "RNA"
  }
  
  if (!ezIsSpecified(param$nreads)) {
    scData$qc.lib <- isOutlier(scData@meta.data[,att_nCounts], log = TRUE, nmads = param$nmad, type = "lower")
  } else {
    scData$qc.lib <- scData@meta.data[,att_nCounts] < param$nreads
  }
  if (!ezIsSpecified(param$ngenes)) {
    scData$qc.nexprs <- isOutlier(scData@meta.data[,att_nGenes], nmads = param$nmad, log = TRUE, type = "lower")
  } else {
    scData$qc.nexprs <- scData@meta.data[,att_nGenes] < param$ngenes
  }
  if (!ezIsSpecified(param$perc_mito)) {
    scData$qc.mito <- isOutlier(scData@meta.data[,"percent_mito"], nmads = param$nmad, type = "higher")
  } else {
    scData$qc.mito <- scData@meta.data[,"percent_mito"] > param$perc_mito
  }
  if (!ezIsSpecified(param$perc_riboprot )) {
    scData$qc.riboprot <- isOutlier(scData@meta.data[,"percent_riboprot"], nmads = param$nmad, type = "higher")
  } else {
    scData$qc.riboprot <- scData@meta.data[,"percent_riboprot"] > as.numeric(param$perc_riboprot)
  }
  
  scData$useCell <- !(scData$qc.lib | scData$qc.nexprs | scData$qc.mito | scData$qc.riboprot)
  
  set.seed(38)
  doubletsInfo <- scDblFinder(GetAssayData(scData, layer="counts")[ , scData$useCell], returnType = "table", clusters=TRUE, BPPARAM = BPPARAM)
  scData$doubletScore <- doubletsInfo[colnames(scData), "score"]
  scData$doubletClass <- doubletsInfo[colnames(scData), "class"]
  scData$qc.doublet <- scData$doubletClass %in% "doublet"
  if (ezIsSpecified(param$keepDoublets) && param$keepDoublets) {
    futile.logger::flog.info("Keeping doublets...")
  } else {
    scData$useCell <- scData$useCell & scData$doubletClass %in% "singlet"
  }
  return(scData)
}

querySignificantClusterAnnotationEnrichR <- function(genesPerCluster, dbs, overlapGeneCutOff = 3, adjPvalueCutOff = 0.001, reportTopN = 5) {
  enrichRout <- list()
  for (cluster in unique(names(genesPerCluster))) {
    enriched <- enrichr(as.character(genesPerCluster[[cluster]]), dbs)
    
    for (db in names(enriched)) {
      enriched_db <- enriched[[db]]
      if (nrow(enriched_db) > 0 && colnames(enriched_db)[1] == "Term"){
        enriched_db$OverlapGenesN <- sub("/.*", "", enriched_db$Overlap) %>% as.numeric()
        enriched_db$Cluster <- cluster
        enriched_db <- enriched_db %>%
          filter(., Adjusted.P.value < adjPvalueCutOff) %>%
          filter(., OverlapGenesN > overlapGeneCutOff) %>%
          head(reportTopN)
        enrichRout[[cluster]][[db]] <- enriched_db[, c("Term", "Cluster", "Overlap", "OverlapGenesN", "Adjusted.P.value", "Odds.Ratio", "Combined.Score")]
      }
    }
  }
  return(enrichRout)
}


computeTFActivityAnalysis <- function(cells, species){
  species <- tolower(species)
  # Retrieve prior knowledge network.
  network <- decoupleR::get_dorothea(organism = species,
                                     levels = c("A", "B", "C"))
  
  # Run weighted means algorithm.
  activities <- decoupleR::run_wmean(mat = as.matrix(GetAssayData(cells)),
                                     network = network,
                                     .source = "source",
                                     .targe = "target",
                                     .mor = "mor",
                                     times = 100,
                                     minsize = 5)
  
  return(activities)
}


computePathwayActivityAnalysis <- function(cells, species){
  species <- tolower(species)
  # Retrieve prior knowledge network.
  network <- decoupleR::get_progeny(organism = species)
  
  # Run weighted means algorithm.
  activities <- decoupleR::run_wmean(mat = as.matrix(GetAssayData(cells)),
                                     network = network,
                                     .source = "source",
                                     .targe = "target",
                                     .mor = "weight",
                                     times = 100,
                                     minsize = 5)
  
  return(activities)
}

##Run Seurat
setwd("/srv/GT/analysis/zajacn/p28443/ScSeurat_NoIntrons_SMRTLinkGTF")
input = ezRead.table("/srv/gstore/projects/p28443/CellRanger_NoIntrons_SMRTLinkGTF/dataset.tsv")
cmDirs <- input$`CountMatrix [Link]`
samples = rownames(input)
args = commandArgs(trailingOnly=TRUE)
print(args)
order = as.integer(args[1])

cmDir=paste0("/srv/GT/analysis/zajacn/", cmDirs[order])
sample=samples[order]

if (!dir.exists(sample)){
  dir.create(sample)
} else {
  print("Dir already exists!")
}

setwd(sample)
param = readRDS(paste0("/srv/gstore/projects/p28443/o30669_ScSeurat_noIntrons_2023-07-26--11-52-54/", sample, "_SCReport/param.rds"))
param$cellbender = FALSE
if (param$cores > 1){
  BPPARAM <- MulticoreParam(workers = param$cores)
} else {
  ## scDblFinder fails with many cells and MulticoreParam
  BPPARAM <- SerialParam() 
}

cts <- Read10X(cmDir, gene.column = 1)
featInfo <- ezRead.table(paste0(cmDir, "/features.tsv.gz"), header = FALSE, row.names = NULL)#, col_names = FALSE)
colnames(featInfo) <- c("gene_id", "gene_name", "type")
featInfo$isMito = grepl( "(?i)^MT-", featInfo$gene_name)
featInfo$isRiboprot = grepl(  "(?i)^RPS|^RPL", featInfo$gene_name)
featInfo$gene_ID = str_replace(featInfo$gene_id, pattern = ".[0-9]+$",replacement = "")
geneAnno <- ezRead.table("/srv/GT/reference/Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_39-2021-12-09/Genes/SMRTLink11.1_annotation/genes_annotation_byGene.txt")
featInfo$isRibosomal <- geneAnno[featInfo$gene_ID, "type"] == "rRNA"
if(any(is.na(featInfo[, "isRibosomal"]))){
  featInfo[, "isRibosomal"][which(is.na(featInfo[, "isRibosomal"]))] <- FALSE
}

if (is.list(cts)){
  cts <- cts$`Gene Expression`
  featInfo <- featInfo[  featInfo$type == "Gene Expression", ]
}
if(param$cellbender){
  rownames(featInfo) <- featInfo$gene_id
  matchingIds <- intersect(rownames(cts), rownames(featInfo))
  cts <- cts[matchingIds,]
  featInfo <- featInfo[matchingIds,]
}

rownames(cts) <- rownames(featInfo) <- gsub("_", "-", uniquifyFeatureNames(ID=featInfo$gene_ID, names=featInfo$gene_name))

scData <- CreateSeuratObject(counts = cts[rowSums2(cts >0) >0, ])
scData@meta.data$Sample <- sample
scData[["RNA"]] <- AddMetaData(object = scData[["RNA"]], metadata = featInfo[rownames(scData), ])
scData$cellBarcode <- sub(".*_", "", colnames(scData))
scData <- addCellQcToSeurat(scData, param=param, BPPARAM = BPPARAM, ribosomalGenes = featInfo[rownames(scData), "isRibosomal"])

## use empty drops to test for ambient
if ("UnfilteredCountMatrix" %in% colnames(input)) {
  rawDir <- input$UnfilteredCountMatrix
  if (file.exists(file.path(rawDir, param$geneCountModel))) {
    rawDir <- file.path(rawDir, param$geneCountModel)
  }
} else {
  # DEPRECATED; all input datasets should specify the path to the unfiltered count matrix
  rawDir <- sub("filtered_", "raw_", cmDir)
}

if (file.exists(rawDir) && rawDir != cmDir){
  if(param$cellbender){
    rawCts <- Read10X_h5(file.path(dirname(cmDir), 'cellbender_raw_seurat.h5'), use.names = FALSE)
  } else {
    rawCts <- Read10X(rawDir, gene.column = 1)
  }
  if (is.list(rawCts)) {
    rawCts <- rawCts$`Gene Expression`
    rawCts <- rawCts[featInfo$gene_id,]
  }
  
  if (is.null(colnames(input)[grepl("SCDataOrigin", colnames(input))])  && 
      is.null(input[rownames(input) == sample,colnames(input)[grepl("SCDataOrigin", colnames(input))]]  == 'BDRhapsody')) {
    rawCts <- rawCts[featInfo$gene_id,]
  }
  
  if(param$cellbender){
    rawCts <- rawCts[featInfo$gene_id,]
  }
  
  if(length(setdiff(rownames(rawCts), featInfo$gene_id)) > 0){
    rawCts <- rawCts[featInfo$gene_id,] 
  }
  stopifnot(rownames(rawCts) == featInfo$gene_id)
  emptyStats <- emptyDrops(rawCts[!featInfo$isMito & !featInfo$isRiboprot, ],
                           BPPARAM=BPPARAM, niters=1e5)
  scData$negLog10CellPValue <- - log10(emptyStats[colnames(scData), "PValue"])
  emptyStats <- emptyDrops(rawCts, BPPARAM=BPPARAM, niters=1e5)
  scData$negLog10CellPValue <- pmin(scData$negLog10CellPValue, -log10(emptyStats[colnames(scData), "PValue"]))
  scData@meta.data$negLog10CellPValue[is.na(scData$negLog10CellPValue)] <- 0
  remove(rawCts)
}
allCellsMeta <- scData@meta.data
scData <- subset(scData, cells=rownames(allCellsMeta)[allCellsMeta$useCell]) # %>% head(n=1000))

## remove low expressed genes
num.cells <- param$cellsFraction * ncol(scData) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
cellsPerGene <- Matrix::rowSums(GetAssayData(scData, layer="counts") >= param$nUMIs)
is.expressed <- cellsPerGene >= num.cells
cellsPerGeneFraction <- data.frame(frac = cellsPerGene/ncol(scData), row.names = rownames(cellsPerGene))
scData <- scData[is.expressed,]

## Add Cell Cycle information to Seurat object as metadata columns
scData <- addCellCycleToSeurat(scData, param$refBuild, BPPARAM)
## Get information on which variables to regress out in scaling/SCT
scData <- seuratStandardSCTPreprocessing(scData, param)
## defaultAssay is now SCT
scData <- seuratStandardWorkflow(scData, param, ident.name="seurat_clusters")

# estimate ambient first
if (ezIsSpecified(param$estimateAmbient) && param$estimateAmbient) {
  scData <- addAmbientEstimateToSeurat(scData, rawDir=rawDir, threads = param$cores)
}

# get markers and annotations
anno <- getSeuratMarkersAndAnnotate(scData, param)

# save markers
markers <- anno$markers

writexl::write_xlsx(markers, path="posMarkers.xlsx")

stopifnot(length(sample) == 1)
clusterInfos <- ezFrame(Sample=sample, Cluster=levels(Idents(scData)), ClusterLabel="")
if (!is.null(anno$singler.results)){
  clusterInfos$SinglerCellType <- anno$singler.results$singler.results.cluster[clusterInfos$Cluster, "pruned.labels"]
}
nTopMarkers <- 10
topMarkers <- markers %>% group_by(cluster) %>%
  slice_max(n = nTopMarkers, order_by = avg_log2FC)
topMarkerString <- sapply(split(topMarkers$gene, topMarkers$cluster), paste, collapse=", ")
clusterInfos[["TopMarkers"]] <- topMarkerString[clusterInfos$Cluster]
clusterInfoFile <- "clusterInfos.xlsx"

writexl::write_xlsx(clusterInfos, path=clusterInfoFile)

output=readRDS(paste0("/srv/gstore/projects/p28443/o30669_ScSeurat_noIntrons_2023-07-26--11-52-54/", sample, "_SCReport/output.rds"))
output$meta = ezFrame('Condition [Factor]' = sample, 
                      Species = "Homo sapiens (human)", 
                      refBuild="Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_39-2021-12-09", 
                      refFeatureFile = "SMRTLink11.1_annotation/genes.gtf", 
                      'Static Report [Link]'=paste0("/srv/GT/analysis/zajacn/p28443/ScSeurat_NoIntrons_SMRTLinkGTF/", sample, "/", "00index.html"), 
                      'SC Cluster Report [File]' = paste0("/srv/GT/analysis/zajacn/p28443/ScSeurat_NoIntrons_SMRTLinkGTF/", sample), 
                      'SC Seurat' = paste0("/srv/GT/analysis/zajacn/p28443/ScSeurat_NoIntrons_SMRTLinkGTF/", sample, "/", "scData.rds"), 
                      'Sample Id [B-Fabric]' = NA, 'Order Id [B-Fabric]' = NA)
cmd <- paste0("cp /home/zajac/Genomics/ezRun/inst/templates/ScSeurat.Rmd .")
ezSystem(cmd)

makeRmdReport(param=param, output=output, scData=scData, allCellsMeta=allCellsMeta, 
              cellsPerGeneFraction = cellsPerGeneFraction, enrichRout=anno$enrichRout, 
              cells.AUC=anno$cells.AUC, singler.results=anno$singler.results, aziResults=anno$aziResults,
              pathwayActivity=anno$pathwayActivity, TFActivity=anno$TFActivity,
              rmdFile = "ScSeurat.Rmd", reportTitle = paste0(param$name, ": ",  sample))