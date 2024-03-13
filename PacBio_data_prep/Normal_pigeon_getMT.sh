cd /srv/GT/analysis/zajacn/p28443/p28443_o31598_1/

#4_D01
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
--log-level INFO \
--log-file pigeon-make-seurat-4_D01.log \
--num-threads 7 \
--dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/45ddb604-f912-4973-8189-ee9ff3deadf3/call-pb_sc_isoseq/pb_sc_isoseq/5ed15e73-b648-42a4-b589-1d60c49ece9b/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
--group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/45ddb604-f912-4973-8189-ee9ff3deadf3/call-pb_sc_isoseq/pb_sc_isoseq/5ed15e73-b648-42a4-b589-1d60c49ece9b/call-pigeon/inputs/2039781561/scisoseq.mapped_transcripts.collapse.group.txt \
--out-dir p28443_o31598_1_4_D01_getMTinfo \
--out-prefix scisoseq \
/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/45ddb604-f912-4973-8189-ee9ff3deadf3/call-pb_sc_isoseq/pb_sc_isoseq/5ed15e73-b648-42a4-b589-1d60c49ece9b/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt

#2_C01
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
--log-level INFO \
--log-file pigeon-make-seurat-2_C01.log \
--num-threads 7 \
--dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/be5bc430-4db5-42d1-b7b0-17e583ccb27c/call-pb_sc_isoseq/pb_sc_isoseq/0f15096a-8931-4c42-af09-b3aabc07c9a3/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
--group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/be5bc430-4db5-42d1-b7b0-17e583ccb27c/call-pb_sc_isoseq/pb_sc_isoseq/0f15096a-8931-4c42-af09-b3aabc07c9a3/call-pigeon/inputs/-1260006956/scisoseq.mapped_transcripts.collapse.group.txt \
--out-dir p28443_o31598_1_2_C01_getMTinfo \
--out-prefix scisoseq \
/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/be5bc430-4db5-42d1-b7b0-17e583ccb27c/call-pb_sc_isoseq/pb_sc_isoseq/0f15096a-8931-4c42-af09-b3aabc07c9a3/call-pbreports_sc_isoseq/inputs/-707225370/scisoseq_classification.filtered_lite_classification.txt

#1_B01
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
--log-level INFO \
--log-file pigeon-make-seurat-1_B01.log
--num-threads 7 \
--dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f66e783a-7493-48ee-adc3-4eb5bccb2eb8/call-pb_sc_isoseq/pb_sc_isoseq/a1aa7b1d-9036-46ea-bb86-43034c06d6e3/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
--group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f66e783a-7493-48ee-adc3-4eb5bccb2eb8/call-pb_sc_isoseq/pb_sc_isoseq/a1aa7b1d-9036-46ea-bb86-43034c06d6e3/call-pigeon/inputs/1193010483/scisoseq.mapped_transcripts.collapse.group.txt \
--out-dir p28443_o31598_1_1_B01_getMTinfo \
--out-prefix scisoseq \
/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f66e783a-7493-48ee-adc3-4eb5bccb2eb8/call-pb_sc_isoseq/pb_sc_isoseq/a1aa7b1d-9036-46ea-bb86-43034c06d6e3/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt

#1_B01_subsampled
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
--log-level INFO \
--log-file pigeon-make-seurat-1_B01_subsampled.log \
--num-threads 7 \
--dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/b890e18b-9b8b-4fbe-a5e0-ec9fe2a3a066/call-pb_sc_isoseq/pb_sc_isoseq/24f9f296-5d83-4557-83ad-d69d8ec5b080/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
--group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/b890e18b-9b8b-4fbe-a5e0-ec9fe2a3a066/call-pb_sc_isoseq/pb_sc_isoseq/24f9f296-5d83-4557-83ad-d69d8ec5b080/call-pigeon/inputs/487465881/scisoseq.mapped_transcripts.collapse.group.txt \
--out-dir 28443_o31598_1_1_B01_getMTinfo_subsampled \
--out-prefix scisoseq \
/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/b890e18b-9b8b-4fbe-a5e0-ec9fe2a3a066/call-pb_sc_isoseq/pb_sc_isoseq/24f9f296-5d83-4557-83ad-d69d8ec5b080/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt

#3cells_combined - rerun all for all types of data
DIR=/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/0682452f-1241-41a2-8774-f67a2004c9ba
DEDUP=$DIR/call-pb_sc_isoseq/pb_sc_isoseq/e3b856b4-e690-4d6c-a024-6a369ed1e8e5/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta
GROUP=$DIR/call-pb_sc_isoseq/pb_sc_isoseq/e3b856b4-e690-4d6c-a024-6a369ed1e8e5/call-pigeon/inputs/-1437824940/scisoseq.mapped_transcripts.collapse.group.txt
FILTERED=$DIR/call-pb_sc_isoseq/pb_sc_isoseq/e3b856b4-e690-4d6c-a024-6a369ed1e8e5/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt
UNFILTERED=$DIR/call-pb_sc_isoseq/pb_sc_isoseq/e3b856b4-e690-4d6c-a024-6a369ed1e8e5/call-pigeon/execution/scisoseq_classification.txt
REASONS=$DIR/call-pb_sc_isoseq/pb_sc_isoseq/e3b856b4-e690-4d6c-a024-6a369ed1e8e5/call-pigeon/execution/scisoseq_classification.filtered_lite_reasons.txt

for i in $(ls -d p28443_o31598_1_threecells*)
do
  if [[ -d "$i" && "$i" =~ unfiltered ]]; then ##for the data not filtered from isoforms
  /srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
    --log-level INFO \
    --log-file pigeon-make-seurat_threecells.log \
    --num-threads 7 \
    --dedup $DEDUP \
    --group $GROUP \
    --out-dir $i \
    --out-prefix scisoseq \
    $UNFILTERED
  elif [[ -d "$i" && "$i" =~ filteredIntraPriming ]]; then ##for the data filtered only from IntraPriming isoforms
  awk 'NR == FNR {x[$1] = 1; next} ($1 in x) {print $0}' <(cat $REASONS | grep -v 'IntraPriming' | awk -F"," '{print $1}') $UNFILTERED > temp_classification.txt ##obtain the isoforms classified as anything else but intrapriming
  cat $FILTERED temp_classification.txt > temp.classification.txt ##add those isoforms to the filtered classification file
  /srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_threecells.log \
  --num-threads 7 \
  --dedup $DEDUP \
  --group $GROUP \
  --out-dir $i \
  --out-prefix scisoseq \
  temp.classification.txt  
  else ##for the filtered data
  /srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \ 
    --log-level INFO \
    --log-file pigeon-make-seurat_threecells.log \
    --num-threads 7 \
    --dedup $DEDUP \
    --group $GROUP \
    --out-dir $i \
    --out-prefix scisoseq \
    $FILTERED; fi
;done
  

