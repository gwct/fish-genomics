# Selection pipeline

From annotated genomes to selection tests

## Intro

Each step of the pipeline is contained within one script, sometimes python sometimes snakemake. For python, some paths may be hard coded. For snakemake, a config file () provides all the paths and will need to be updated for your file system. There is also a snakemake profile () that allocates resources for your given cluster.

## 0. Starting point

Unaligned orthologous coding sequences (codons) from each genome. While your ortho predictions are probably done with peptide sequences, you will need to retrieve the corresponding CDS sequences to use this pipeline.

## 1. Align orthologous CDS sequences with Guidance

