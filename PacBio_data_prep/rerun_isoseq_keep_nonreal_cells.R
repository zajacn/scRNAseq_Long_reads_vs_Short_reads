sampleInfo <- ezRead.table(file="2024-05-22-sampleInfo.tsv")

ezSystem("cd /scratch/zajacn/p28443/all_cells/")

for (i in unique(sampleInfo$Name)){
  df = sampleInfo[sampleInfo$Name ==i,]
  sampleID = str_replace(df$PbTubeId, "/", "_")
  dir = str_remove(sampleID,"o")
  pacbiodir = sapply(str_split(df$IsoSeqBam, "/"), .subset, 7)
  smrtdir="/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin"
  
  #ezSystem(paste("cd ", dir))
  set_of_reads=paste0("/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/",pacbiodir, "/",
                      list.files(
                        file.path("/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/",pacbiodir), 
                        recursive = TRUE, pattern = "scisoseq.5p--3p.tagged.refined.corrected.sorted.transcriptset.xml")[3])
  
  print(set_of_reads)
  print(file.exists(set_of_reads))
  
  #Deduplicate and keep real cells
  cmd1 = paste(file.path(smrtdir,"isoseq groupdedup"), 
               "--log-level INFO", 
               "--log-file isoseq_dedup.log", 
               "--keep-non-real-cells", 
               " --alarms alarms.json", "-j 7", 
               "`readlink -f", set_of_reads,"`", 
               "scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.bam")
  print(cmd1)
  #ezSystem(cmd1)
  
  cmd2 = paste(file.path(smrtdir,"dataset"),
               "--strict" ,"--log-level INFO" ,
               "--log-file dataset.log",
               "create",
               "--generateIndices", 
               "--type TranscriptSet" ,
               "scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.transcriptset.xml" ,
               "scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.bam")
  print(cmd2)
  #ezSystem(cmd2)
  
  ezSystem(paste("echo -n 'scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.bam' > unmapped_bam.fofn"))
  ezSystem(paste("echo -n 'scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta' > dedup_fasta.fofn"))
  
  reference=paste0("/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/", 
                   pacbiodir, "/", 
                   list.files(
                     file.path("/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/", 
                               pacbiodir), recursive = TRUE, pattern = "Human_hg38_Gencode_v39.referenceset.xml")[1]
                   )
  print(reference)
  print(file.exists(reference))
  #Map reads to reference
  cmd3 = paste(file.path(smrtdir,"pbmm2 align"), 
               "--log-level INFO",
               "--log-file pbmm2.log",
               "-j 7" ,
               "--alarms alarms.json" ,
               "--preset ISOSEQ" ,"--sort" ,
               "--report-json mapping_stats.report.json" ,
               "`readlink -f", reference,"`",
               "scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.transcriptset.xml",
               "scisoseq.mapped.transcriptalignmentset.xml")
  #ezSystem(cmd3)
  print(cmd3)
  ezSystem(paste("echo -n 'scisoseq.mapped.bam' > mapped_bam.fofn"))
  ezSystem(paste("echo -n 'scisoseq.mapped.bam.bai' > mapped_bam_bai.fofn"))
  
  
  #Collapse transcripts - using criteria from SMRTLink v11.1 and keeping non real cells
  cmd4 = paste(file.path(smrtdir,"isoseq collapse"),
               "--log-level INFO", 
               "--log-file isoseq_collapse.log", 
               "-j 7",
               "--keep-non-real-cells", 
               "--max-5p-diff 1000", 
               "--max-3p-diff 100", 
               "--alarms alarms.json", 
               "scisoseq.mapped.transcriptalignmentset.xml", 
               "scisoseq.mapped_transcripts.collapse.gff")
  print(cmd4)
  #ezSystem(cmd4)
  
  cmd5 = paste(file.path(smrtdir,"dataset"),
               "--log-level INFO", 
               "--log-file dataset.log", 
               "absolutize", 
               "--output-xml absolutized.referenceset.xml", 
               reference)
  #ezSystem(cmd5)
  print(cmd5)
  #Start the pigeon workflow
  cmd6 = paste(file.path(smrtdir,"pigeon sort"), 
               "--log-level INFO",
               "--log-file pigeon-sort.log",
               "-o scisoseq_transcripts.sorted.gff",
               "scisoseq.mapped_transcripts.collapse.gff")
  #ezSystem(cmd6)
  print(cmd6)

  cmd7 = paste(file.path(smrtdir,"pigeon classify"), 
               "--log-level INFO", 
               "--log-file pigeon-classify.log", 
               "-j 7", 
               "--out-dir .", 
               "--out-prefix scisoseq", 
               "--flnc scisoseq.mapped_transcripts.collapse.abundance.txt", 
               "--ref absolutized.referenceset.xml", 
               "scisoseq_transcripts.sorted.gff")
  #ezSystem(cmd7)
  print(cmd7)
  cmd8 = paste(file.path(smrtdir,"pigeon report"), 
               "--log-level INFO" ,
               "--log-file pigeon-report.log" ,
               "scisoseq_classification.txt" ,
               "scisoseq_saturation.txt")
  #ezSystem(cmd8)
  print(cmd8)
  cmd9 = paste(file.path(smrtdir,"pigeon make-seurat"),
               "--keep-ribo-mito-genes",
               "--log-level INFO" ,
               "--log-file pigeon-make-seurat.log" ,
               "--num-threads 7",
               "--annotations absolutized.referenceset.xml" ,
               "--dedup scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta" ,
               "--group scisoseq.mapped_transcripts.collapse.group.txt" ,
               "--out-dir o31598_5_pigeon",
               "--out-prefix scisoseq" ,
               "scisoseq_classification.txt")
  print(cmd9)
  #ezSystem(cmd9)
  #cd("/scratch/zajacn/p28443/all_cells")
}