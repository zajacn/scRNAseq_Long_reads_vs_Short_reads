#Download blaze
apptainer pull docker://quay.io/biocontainers/blaze2:2.5.1--pyhdfd78af_0

#Convert segmented reads to fastq
/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/bam2fastq segmented.bam -c 9 -o $i.segmented.fastq.gz
#For cells 1 and 2 combine all 3 segmented.bam files

#Run blaze - expected cells = Number from Illumina on fgcz-c-052
apptainer run --bind /scratch/zajac/p28443 blaze2_2.5.1--pyhdfd78af_0.sif blaze --expect-cells 588 --output-prefix 31598_3 --threads 7 ./31598_3/
apptainer run --bind /scratch/zajac/p28443 blaze2_2.5.1--pyhdfd78af_0.sif blaze --expect-cells 1664  --output-prefix 31598_4 --threads 7 ./31598_4/
apptainer run --bind /scratch/zajac/p28443 blaze2_2.5.1--pyhdfd78af_0.sif blaze --expect-cells 1155  --output-prefix 31598_5 --threads 7 ./31598_5/
apptainer run --bind /scratch/zajac/p28443 blaze2_2.5.1--pyhdfd78af_0.sif blaze --expect-cells 827  --output-prefix 31598_1 --threads 7 ./31598_1/
apptainer run --bind /scratch/zajac/p28443 blaze2_2.5.1--pyhdfd78af_0.sif blaze --expect-cells 549  --output-prefix 31598_2 --threads 7 ./31598_2/
