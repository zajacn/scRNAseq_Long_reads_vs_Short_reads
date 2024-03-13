cd /srv/GT/analysis/zajacn/p28443/p28443_o31598_5/

#4_E01
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_4_E01.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-pigeon/inputs/845432975/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir  p28443_o31598_5_4_E01_getMTinfo \
  --out-prefix scisoseq \
  /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt

#4_E01_unfiltered
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_4_E01.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-pigeon/inputs/845432975/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir  p28443_o31598_5_4_E01_getMTinfo_unfiltered \
  --out-prefix scisoseq \
  /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-pigeon/execution/scisoseq_classification.txt

#4_E01_filteredIntraPriming
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_4_E01.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-pigeon/inputs/845432975/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir  p28443_o31598_5_4_E01_getMTinfo_filteredIntraPriming \
  --out-prefix scisoseq \
  temp.classification.txt



awk 'NR == FNR {x[$1] = 1; next} ($1 in x) {print $0}' <(cat /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-pigeon/execution/scisoseq_classification.filtered_lite_reasons.txt | grep -v 'IntraPriming' | awk -F"," '{print $1}') /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-pigeon/execution/scisoseq_classification.txt > temp_classification.txt
cat /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/e3fea0f1-c152-4ffc-862f-d0e2b5af03fc/call-pb_sc_isoseq/pb_sc_isoseq/1c5cb782-8770-4e00-b4f2-a5523b3f8258/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt temp_classification.txt > temp.classification.txt
