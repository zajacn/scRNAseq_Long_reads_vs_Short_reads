cd /srv/GT/analysis/zajacn/p28443/p28443_o31598_4/

#3_D01
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_3_D01.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-pigeon/inputs/-2129936381/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir p28443_o31598_4_3_D01_getMTinfo \
  --out-prefix scisoseq \
  /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt

#3_D01_unfiltered
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_3_D01.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-pigeon/inputs/-2129936381/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir p28443_o31598_4_3_D01_getMTinfo_unfiltered \
  --out-prefix scisoseq \
  /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-pigeon/execution/scisoseq_classification.txt

#3_D01_filteredIntraPriming
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_3_D01.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-pigeon/inputs/-2129936381/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir p28443_o31598_4_3_D01_getMTinfo_filteredIntraPriming \
  --out-prefix scisoseq \
  temp.classification.txt


awk 'NR == FNR {x[$1] = 1; next} ($1 in x) {print $0}' <(cat /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-pigeon/execution/scisoseq_classification.filtered_lite_reasons.txt | grep -v 'IntraPriming' | awk -F"," '{print $1}') /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-pigeon/execution/scisoseq_classification.txt > temp_classification.txt
cat /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/9444586f-5bc1-47f3-be9b-1e80291de67a/call-pb_sc_isoseq/pb_sc_isoseq/256183fd-b4f7-41da-b225-4ef0cbac0c7e/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt temp_classification.txt > temp.classification.txt
