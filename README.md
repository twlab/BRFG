# Benchmarking references for functional genomics
This repo contains the code used in performing the study.

## Compute 
Two compute systems were used in this study.
The WashU HTCF and WashU Compute1 systems. HTCF is a slurm based system. Compute1 is a docker based system using LSF for job submission.

- Workflows
  - ATAC-seq: [AIAP](AIAP)
  - RNA-seq: ENCODE
  - Hi-C: ENCODE
  - WGBS:
- Post-processing analysis
  - Read tracing: [Read_Tracing](Read_Tracing)
  - QC analysis:
    - [Plot_Metrics.Rmd](Analysis_R_Scripts/HPRC_ATACseq/Plot_Metrics.Rmd)
    - [PlotMetrics.R](Analysis_R_Scripts/HPRC_RNAseq/PlotMetrics.R)
    - [Plot_Metrics_and_Outputs.Rmd](Analysis_R_Scripts/HPRC_HiC/Plot_Metrics_and_Outputs.Rmd)
    - [PlotMetrics.R](Analysis_R_Scripts/HPRC_WGBS/PlotMetrics.R)
  - Functional estimate analysis
    - [ComparePeaks.Rmd](Analysis_R_Scripts/HPRC_ATACseq/ComparePeaks.Rmd)
    - [CompareExpressionByStrategy.Rmd](Analysis_R_Scripts/HPRC_RNAseq/CompareExpressionByStrategy.Rmd)
    - [Compare_LinearOnly_Latest.R](Analysis_R_Scripts/HPRC_WGBS/Compare_LinearOnly_Latest.R)
- Variant enrichment analysis
  - [HPRC_DV_Enrichment](Analysis_R_Scripts/HPRC_DV_Enrichment)
  - [HPRC_SV_Enrichment](Analysis_R_Scripts/HPRC_SV_Enrichment)
  - [Pangenome_structure_diversity](Analysis_R_Scripts/Pangenome_structure_diversity)
- Overrepresentation analysis
  - [ORA](Analysis_R_Scripts/ORA)
- Liftover
  - [Liftover](Liftover)