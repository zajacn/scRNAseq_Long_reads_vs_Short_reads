# scRNAseq_Long_reads_vs_Short_reads
This repository contains a code for a publication entitled: Comparison of Single-cell Long-read and Short-read Transcriptome Sequencing of Patient-derived Organoid Cells of ccRCC: Quality Evaluation of the MAS-ISO-seq Approach: https://doi.org/10.1093/nargab/lqaf089

bioRxiv_Scripts refer to the bioRxiv version of this article: https://doi.org/10.1101/2024.03.14.584953.


The repository contains the following set of scripts/reports:

- Illumina Data preparation: Data has been processed with CellRanger v7.2.    

- PacBio Data preparation: Iso-seq and Pigeon have been run as part of the SMRTLink v11.1 software available open source: see https://www.pacb.com/support/software-downloads/.     
  - In order to recover cell barcodes that were not called real cells we used the script: rerun_isoseq_keep_nonreal_cells.R.     
  - In order to recover transcripts and reads mapping to mitochondrial reads we reran Pigeon using the *_pigeon_getMT.sh scripts.     

- Whole data summary computed from mapped data.
  - Cell statistics (mitochondrial content, TSO and UMI counts) computed from PacBio and Illumina bam files: compileMITOstats.R, compileTSOPAInfo.R, compileUMIcounts.R

- Cell calling: 
  - Empty Droplet analysis is not part of the Iso-seq workflow, we ran BLAZE (see https://github.com/shimlab/BLAZE) for that: blaze.sh.     
  - Analysis of cell barcodes called real cells uniquely in PacBio by Iso-seq: 2024-11-05-cells_unique_to_PacBio_anaylsis.Rmd.    
  - Analysis of cell barcodes called real cells uniquely in Illumina by CellRanger (low RNA content cells):2024-12-15-lowRNAcontentcells_PacBiovsIllumina.Rmd.   

- Comparison of a subset of mapped data between Illumina and PacBio:
  - Convert bam files to dataframes summarising the information on common and unique cell barcode- UMI tags: compileMatchedReads.R.
  - Reports/Figures for the data subsampled 3 times (for cell barcodes starting with TG, GC and AA): 2024-12-09-summarise-bam-subset-TG.Rmd, 2024-12-09-summarise-bam-subset-GC.Rmd, 2024-12-09-summarise-bam-subset.Rmd.
  - A table summarising the differences in wet lab and data analysis steps between the two technologies: Methods_Comparison_table.txt.

- Comparison of the assembled gene count matrix from Illumina data using CellRanger workflow and from PacBio data using Iso-seq + pigeon workflow:
  - Report/Figures: 2025-01-13-gene-count-matrix-comparison-pacbiounfiltered-Illumina.Rmd.

- Comparison of PacBio data filtered and unfiltered of artefactual isoforms as part of the Iso-seq + Pigeon workflow (RT-switching, LCNC, IntraPriming):
  - Report/Figures: 2025-01-19-PacBio-filtering-analysis.Rmd.

Note: All scripts contain hardcoded paths to the location of each dataset on our system therefore they would need to be all modified for running them on a different location. 





