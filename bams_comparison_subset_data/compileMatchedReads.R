#!/usr/bin/Rscript
options(error=recover)
library(ezRun)
library(tidyverse)
library(parallel)
library(dplyr)
library(Rsamtools)
library(GenomicAlignments)

sample_file = read_tsv("../2024-05-22-sampleInfo.tsv")
sample_file = sample_file[5,]
workDir <- "/srv/GT/analysis/zajacn/p28443/"

compileMatchedReads_20240603 <- function(sampleInfo= sample_file, cbStart="AA", fastOverlapOnly=FALSE, nCores=4,
                                         scratchDir="/srv/GT/analysis/zajacn/p28443/subset_bams_072024"){
  
  if (!file.exists(scratchDir)){
    dir.create(scratchDir)
  }
  i <- 1
  sampleInfo$pbGaFile <- paste0(scratchDir, "/", sampleInfo$Name, "-", cbStart, "-pbGa.qsd")
  sampleInfo$illGaFile <- paste0(scratchDir, "/", sampleInfo$Name, "-", cbStart, "-illGa.qsd")
  sampleInfo$illUnmappedFile <- paste0(scratchDir, "/", sampleInfo$Name, "-", cbStart, "-illUnmapped.qsd")
  for (i in seq_along(sampleInfo$Name)){
    smi <- sampleInfo[i, ]
    message(smi$Name)
    if (!file.exists(smi$pbGaFile)){
      tmpSam <- paste0(scratchDir, "/", smi$Name, "-", cbStart, ".sam")
      illDir <- paste0("/srv/gstore/projects/", smi$CellRangerNoIntron)
      illBam <- sub(".sam$", ".bam", tmpSam)
      if (!file.exists(illBam)){
        ezSystem(paste0("/usr/local/ngseq/packages/Tools/samtools/1.11/bin/samtools view -H ", 
                        illDir, "/possorted_genome_bam.bam > ",  tmpSam))
        ezSystem(paste0("/usr/local/ngseq/packages/Tools/samtools/1.11/bin/samtools view ", 
                        illDir, "/possorted_genome_bam.bam| grep CB:Z:", cbStart, " >>",  tmpSam))
        ezSystem(paste("/usr/local/ngseq/packages/Tools/samtools/1.11/bin/samtools view -b", tmpSam, ">", illBam))
        ezSystem(paste("/usr/local/ngseq/packages/Tools/samtools/1.11/bin/samtools index", illBam))
        file.remove(tmpSam)
      }
      ## read the illumina alignments
      sbp <- ScanBamParam(what =  c("qname", "seq"), tag = c("GN", "xf", "CB", "UB", "UR", "ts", "pa", "RE", "NH"))
      illGa <- readGAlignments(illBam, param = sbp)
      mcols(illGa)$CB <- sub("-1", "", mcols(illGa)$CB)
      ## replace the undefined corrected barcodes with the raw barcode
      toReplace <- is.na(mcols(illGa)$UB)
      mcols(illGa)$UB[toReplace] <- mcols(illGa)$UR[toReplace]
      mcols(illGa)$UR <- NULL
      
      sbp <- ScanBamParam(what = c("qname", "seq"), tag = c("xf", "CB", "UB", "ts", "pa"), 
                          scanBamFlag(isUnmappedQuery = TRUE))
      illBamList <- scanBam(illBam, param = sbp)[[1]]
      illBamList$tag$CB <- sub("-1", "", illBamList$tag$CB)
      illUnmapped <- cbind(qname=illBamList$qname, data.frame(illBamList$tag))
      illUnmapped$gc <- letterFrequency(illBamList$seq, letters="GC", as.prob = FALSE)[ ,"G|C"] /width(illBamList$seq)
      illUnmapped$qlength <- width(illBamList$seq)
      
      #Load the Pacbio data and keep only same cells in all datasets
      sbp <- ScanBamParam(what = c("qname", "seq"), tag = c( "CB", "XM"), tagFilter=list(rc=1))
      pbGa <- readGAlignments(smi$IsoSeqBam, param = sbp)
      ## reverse complement the barcode to match illumina
      mcols(pbGa)$CB <- mcols(pbGa)$CB %>% DNAStringSet() %>% reverseComplement() %>% as.character()
      
      ## intersect the alignments with the illumina alignments by cell barcode; pb alignments on contain pb real cells
      pbGa <- pbGa[mcols(pbGa)$CB %in% mcols(illGa)$CB]
      illGa <- illGa[mcols(illGa)$CB %in% mcols(pbGa)$CB]
      illUnmapped <- illUnmapped[illUnmapped$CB %in% mcols(pbGa)$CB, ]
      mcols(pbGa)$tagId <- paste(mcols(pbGa)$CB,
                                 as.character(reverseComplement(DNAStringSet(mcols(pbGa)$XM))))
      mcols(illGa)$tagId <- paste(mcols(illGa)$CB, mcols(illGa)$UB)
      illUnmapped$tagId <- paste(illUnmapped$CB, illUnmapped$UB)
      
      pbAnno <- data.table::fread(file.path(smi$PigeonUnfiltered, "scisoseq.annotated.info.csv"), sep="\t")
      pbAnno <- pbAnno[pbAnno$BCrev %in% mcols(illGa)$CB, ]
      
      pbReasons <- data.table::fread(smi$isoSeqFiltReasons, sep=",", skip="filtered_isoform,filter")
      pbAnno$quality <- pbReasons[match(pbAnno$pbid, pbReasons$filtered_isoform), "filter"] %>% unlist() %>% replace_na("Pass")
      stopifnot(pbAnno$id %in% mcols(pbGa)$qname)
      idx <- match(mcols(pbGa)$qname, pbAnno$id)
      mcols(pbGa)$category <- pbAnno$category[idx] %>% replace_na("unannotated")
      mcols(pbGa)$gene <- pbAnno$gene[idx]
      mcols(pbGa)$transcriptLength <- pbAnno$length[idx]
      mcols(pbGa)$quality <- pbAnno$quality[idx] %>% replace_na("Unknown")
      
      commonTagIds <- intersect(mcols(pbGa)$tagId, mcols(illGa)$tagId)
      mcols(illGa)$pbStatus <- ifelse(mcols(illGa)$tagId %in% commonTagIds, "pbDetected", "pbAbsent")
      mcols(pbGa)$illStatus <- ifelse(mcols(pbGa)$tagId %in% commonTagIds, "illDetected", "illAbsent")
      mcols(pbGa)$illUnmapped <- ifelse(mcols(pbGa)$tagId %in% illUnmapped$tagId, "illUnmapped", NA)
      mcols(pbGa)$gc <- letterFrequency(mcols(pbGa)$seq, letters="GC", as.prob = FALSE)[ ,"G|C"] /width(mcols(pbGa)$seq)
      mcols(pbGa)$qlength <- width(mcols(pbGa)$seq)
      mcols(illGa)$gc <- letterFrequency(mcols(illGa)$seq, letters="GC", as.prob = FALSE)[ ,"G|C"] /width(mcols(illGa)$seq)
      mcols(illGa)$qlength <- width(mcols(illGa)$seq)
      
      
      mcols(pbGa)$seq <- NULL
      mcols(illGa)$seq <- NULL
      ## subset for faster analysis:
      ### commonTagIds <- head(commonTagIds, 1000)
      if (!fastOverlapOnly){
        ## turn int GRangesList
        job <- ezJobStart("map")
        pbRgList <- as(pbGa[mcols(pbGa)$tagId %in% commonTagIds], "GRangesList")
        ezWriteElapsed(job, status="Grangeslist")
        # ## regroup by tagid
        pbTagId <- mcols(pbRgList)$tagId
        mcols(pbRgList) <- NULL
        pbRgs <- unlist(pbRgList)
        pbRgs$tagId <- rep(pbTagId, elementNROWS(pbRgList))
        pbRgList <- split(pbRgs, pbRgs$tagId) ## this is a GRangesList
        ## group pbAlignments by tagId and discard the mcols; this is a normal array/list -- insanely slow
        # pbRgList <- as(pbGa, "GRangesList")
        # mcols(pbRgList) <- NULL
        # pbRgList <- tapply(pbRgList[1:10000], mcols(pbGa)$tagId[1:10000], unlist, simplify = FALSE)
        ezWriteElapsed(job, status="Grangeslist regrouped")
        useIlluGa <- mcols(illGa)$tagId %in% commonTagIds
        myIllTagId <- mcols(illGa)$tagId[useIlluGa]
        illCommRgList <-  as(illGa[useIlluGa], "GRangesList") # [mcols(illGa)$tagId %in% names(pbRgList)]
        mcols(illCommRgList) <- NULL
        ezWriteElapsed(job, status=paste("illu Grangeslist subsetted to", length(illCommRgList)))
        
        ### computing the overlap is the slow part
        ### using mapply would be slower than mclapply
        # hasOverlap <- mapply(function(x,y){any(overlapsAny(x,y, minoverlap=5, type="any"))},
        #               illCommRgList, pbRgList[myIllTagId])
        ### also simplifying the illumina ranges with: endoapply(illCommRgList, GenomicRanges::reduce)
        ### is not faster
        hasOverlap <- mclapply(1:length(illCommRgList), function(i){ 
          if (i %% 100000 == 0) ezWriteElapsed(job, status=paste("mapped", i))
          any(overlapsAny(illCommRgList[[i]], pbRgList[[myIllTagId[i]]], minoverlap=5, type="any"))
        }, mc.cores=nCores)
        hasOverlap <- unlist(hasOverlap)
        ezWriteElapsed(job, status="mapping done")
        mcols(illGa)$overlapsPbAlign <- NA
        mcols(illGa)$overlapsPbAlign[useIlluGa] <- ifelse(hasOverlap, "hasOverlapPb", "noOverlapPb")
        
        gc()
        ezWriteElapsed(job, status="cleaned up")
      }
      # fast first version which ignores the cigar information; and maps to the first Pb alignment of the tagId
      job <- ezJobStart("fast map")
      useIlluGa <- mcols(illGa)$tagId %in% commonTagIds
      illCommRgs <- granges(illGa[useIlluGa], use.mcols=FALSE)
      pbCommRgs <- granges(pbGa, use.mcols=FALSE)[match(mcols(illGa)$tagId[useIlluGa], mcols(pbGa)$tagId)]
      ovlStatus <- poverlaps(illCommRgs, pbCommRgs,
                             type="any", minoverlap=10) %>% as.vector()
      mcols(illGa)$overlapsPbAlignRange <- NA
      mcols(illGa)$overlapsPbAlignRange[useIlluGa] <- ifelse(as.vector(ovlStatus), "hasOverlapPb", "noOverlapPb")
      ezWriteElapsed(job, status="done")
      # as.vector(ovlStatus[match(mcols(illCommRgs)$tagId, names(ovlStatus)])
      # illCommRgs <- granges(illGa[match(commonTagIds, mcols(illGa)$tagId)], use.mcols=TRUE)
      # pbCommRgs <- granges(pbGa[match(commonTagIds, mcols(pbGa)$tagId)], use.mcols=TRUE)
      # # 
      # ovlStatus <- poverlaps(illCommRgs, pbCommRgs,
      #                        type="any", minoverlap=10) %>% as.vector()
      # # ### summarize duplicates by tagId
      # ovlStatus <- tapply(ovlStatus, mcols(illCommRgs)$qname, any)
      # mcols(illGa)$overlapsPbAlignRange <- as.vector(ovlStatus[mcols(illGa)$qname])
      
      isMultiMapper <- mcols(pbGa)$qname %in% mcols(pbGa)$qname[duplicated(mcols(pbGa)$qname)]
      mcols(pbGa)$multiMapping <-  ifelse(isMultiMapper, "Multi mapping", "Uniquely mapping")
      mcols(illGa)$umiCount = ifelse(mcols(illGa)$xf == "25", "Counted", "NotCounted")
      mcols(illGa)$Sample <- smi$Name
      mcols(pbGa)$Sample <- smi$Name
      qs::qsave(illGa, file=smi$illGaFile)
      qs::qsave(pbGa, file=smi$pbGaFile)
      qs::qsave(illUnmapped, file=smi$illUnmappedFile)
      remove(illGa, pbGa, illCommRgs, pbCommRgs)
      gc()
    }
  }
  return(sampleInfo)
}

sampleInfo <- compileMatchedReads_20240603(sampleInfo = sample_file, 
                                                 scratchDir="/srv/GT/analysis/zajacn/p28443/subset_bams_072024")
