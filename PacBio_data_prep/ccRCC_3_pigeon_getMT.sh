cd /srv/GT/analysis/zajacn/p28443/p28443_o31598_3/

#2_C01
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_2_C01.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-pigeon/inputs/-170257495/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir p28443_o31598_3_2_C01_getMTinfo \
  --out-prefix scisoseq \
  /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt

#2_C01_unfiltered
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_2_C01.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-pigeon/inputs/-170257495/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir p28443_o31598_3_2_C01_getMTinfo_unfiltered \
  --out-prefix scisoseq \
  /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-pigeon/execution/scisoseq_classification.txt


#2_C01_filteredIntraPriming
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_2_C01.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-pigeon/inputs/-170257495/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir p28443_o31598_3_2_C01_getMTinfo_filteredIntraPriming \
  --out-prefix scisoseq \
  temp.classification.txt


awk 'NR == FNR {x[$1] = 1; next} ($1 in x) {print $0}' <(cat /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-pigeon/execution/scisoseq_classification.filtered_lite_reasons.txt | grep -v 'IntraPriming' | awk -F"," '{print $1}') /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-pigeon/execution/scisoseq_classification.txt > temp_classification.txt
cat /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9e2e5e5d-ab85-4f99-8f8e-05f43b2057ca/call-pb_sc_isoseq/pb_sc_isoseq/2f4e813e-e9ba-4ace-a550-2f931e57035d/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt temp_classification.txt > temp.classification.txt
