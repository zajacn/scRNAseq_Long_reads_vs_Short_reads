cd /srv/GT/analysis/zajacn/p28443/p28443_o31598_2/

#Run pigeon on each SMRTcell seaprately - check for batch effects

#3_D01
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
--log-level INFO \
--log-file pigeon-make-seurat-3_D01.log \
--num-threads 7 \
--dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/3d87fb56-4d75-47be-98c4-91ff3b9d97a9/call-pb_sc_isoseq/pb_sc_isoseq/d1715c56-c261-4355-8641-6c713d4d8a87/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
--group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/3d87fb56-4d75-47be-98c4-91ff3b9d97a9/call-pb_sc_isoseq/pb_sc_isoseq/d1715c56-c261-4355-8641-6c713d4d8a87/call-pigeon/inputs/-1677665194/scisoseq.mapped_transcripts.collapse.group.txt \
--out-dir p28443_o31598_1_3_D01_getMTinfo \
--out-prefix scisoseq \
/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/3d87fb56-4d75-47be-98c4-91ff3b9d97a9/call-pb_sc_isoseq/pb_sc_isoseq/d1715c56-c261-4355-8641-6c713d4d8a87/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt

#3_C01
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
--log-level INFO \
--log-file pigeon-make-seurat-3_C01.log \
--num-threads 7 \
--dedup  /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/fd5d2bd6-7e2b-4076-afff-d6206865f258/call-pb_sc_isoseq/pb_sc_isoseq/0f14f86b-ad8d-418f-afb3-da76adc095a5/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
--group  /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/fd5d2bd6-7e2b-4076-afff-d6206865f258/call-pb_sc_isoseq/pb_sc_isoseq/0f14f86b-ad8d-418f-afb3-da76adc095a5/call-pigeon/inputs/-1082556109/scisoseq.mapped_transcripts.collapse.group.txt \
--out-dir p28443_o31598_1_3_C01_getMTinfo \
--out-prefix scisoseq \
/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/fd5d2bd6-7e2b-4076-afff-d6206865f258/call-pb_sc_isoseq/pb_sc_isoseq/0f14f86b-ad8d-418f-afb3-da76adc095a5/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt

#1_B01
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
--log-level INFO \
--log-file pigeon-make-seurat-1_B01.log \
--num-threads 7 \
--dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/54642992-8a09-448b-9424-5129adeb3959/call-pb_sc_isoseq/pb_sc_isoseq/76f446f3-8bfa-49af-8683-84e1e0236f21/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
--group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/54642992-8a09-448b-9424-5129adeb3959/call-pb_sc_isoseq/pb_sc_isoseq/76f446f3-8bfa-49af-8683-84e1e0236f21/call-pigeon/inputs/-728606909/scisoseq.mapped_transcripts.collapse.group.txt  \
--out-dir p28443_o31598_2_1_B01_getMTinfo_unfiltered \
--out-prefix scisoseq \
/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/54642992-8a09-448b-9424-5129adeb3959/call-pb_sc_isoseq/pb_sc_isoseq/76f446f3-8bfa-49af-8683-84e1e0236f21/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt

#Run pigeon on all cells together for final analysis

#3cells_combined
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_threecells.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-pigeon/inputs/2072308913/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir p28443_o31598_2_threecells_getMTinfo \
  --out-prefix scisoseq \
  /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt

#3cells_combined_unfiltered
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_threecells.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-pigeon/inputs/2072308913/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir p28443_o31598_2_threecells_getMTinfo_unfiltered \
  --out-prefix scisoseq \
  /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-pigeon/execution/scisoseq_classification.txt


#3cells_combined_filteredIntraPriming
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat --keep-ribo-mito-genes \
  --log-level INFO \
  --log-file pigeon-make-seurat_threecells.log \
  --num-threads 7 \
  --dedup /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
  --group /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-pigeon/inputs/2072308913/scisoseq.mapped_transcripts.collapse.group.txt \
  --out-dir p28443_o31598_2_threecells_getMTinfo_filteredIntraPriming \
  --out-prefix scisoseq \
  temp.classification.txt

awk 'NR == FNR {x[$1] = 1; next} ($1 in x) {print $0}' <(cat /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-pigeon/execution/scisoseq_classification.filtered_lite_reasons.txt | grep -v 'IntraPriming' | awk -F"," '{print $1}') /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-pigeon/execution/scisoseq_classification.txt > temp_classification.txt
cat /misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/pb_segment_reads_and_sc_isoseq/f9489590-ac62-476a-b05f-3ecc2da7ef32/call-pb_sc_isoseq/pb_sc_isoseq/bde6615c-fab1-4d42-b17e-62a8470bc85b/call-pigeon/execution/scisoseq_classification.filtered_lite_classification.txt temp_classification.txt > temp.classification.txt
