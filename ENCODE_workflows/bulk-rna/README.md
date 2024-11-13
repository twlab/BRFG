# Bulk RNA Benchmarking

Files, scripts, and notes for processing long HiC.

## Docker

A custom docker was built from the _Dockerfile_ and uploaded to __ebelter/mgi:bulk-rna__. This docker was generally taken from the ENCODE Dockerfile, but modified per below:

* samtools version 1.16
* STAR version 2.7.9a
* use args to store versions for simplicity

## Build Indexes Pipeline

This pipeline combines the ENCODE RNA build_genome_index and merge_anno pipelines into one workflow. Additionally, it creates a fasta of transcritpts. Some of the _scripts_ were modified to make tRNA annotations optional, and to utilize advancements in STAR v2.9.1.

### Inputs

See _build-idxs.inputs.json_ for workflow inputs.

### Steps

See _build-idxs.steps_ for task calls with inputs and outputs.

### Outputs

See _build-idxs.outputs.yaml_ for outputs that are copied from the pipeline run.

## ENCODE Bulk RNA Pipeline

1.2.4 https://github.com/ENCODE-DCC/rna-seq-pipeline

WDL was modified to ... FIXME

### Inputs

See _rna-seq-pipeline.inputs.json_ for workflow inputs.

### Steps

See _rna-seq-pipeline.steps_ for task calls with inputs and outputs.

### Outputs

See _rna-seq-pipeline.outputs.yaml_ for outputs that are copied from the pipeline run.

