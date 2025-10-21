
This repository contains code used to process, analyze and visualise Plasmodium falciparum genomic data from raw fastq files as described in the manuscript [Convergent Evolution of Artemisinin and Chloroquine Resistance in Ethiopian Plasmodium falciparum Parasites](https://verixiv.org/articles/2-162)

Raw sequencing data used in this analysis can be downloaded from [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJEB87415)

Bash scripts scripts rely on functioning installations of bwa, python, samtools, gatk, bcftools, and snpEff.<br />
The analysis was done with the following versions: bwa version 0.7.18-r1243, python version 3.12.7, samtools version 1.22.1, gatk v4.1.4.1, bcftools version 1.22, snpEff version 5.3a. <br />
Packages can be installed using: `conda create -n env bwa=0.7.18 samtools=1.22 bcftools=1.22 snpeff=5.3a gatk4=4.1.4.1 python=3.12` (10-20 minutes typical installation time)

R scripts in R version 4.4.2 and using the following packages: ape (5.8.1), pegas (1.3), msa (1.38.0), dplyr (1.1.4), psych (2.4.5), officer (0.6.7), flextable (0.9.7), tidyr (1.3.1), readxl (1.4.5), purrr (1.0.4), writexl (1.5.1), reshape2 (1.4.4), tidyverse (2.0.0), ggplot2 (3.5.2), scales (1.4.0), sf (1.0.20), grid (4.4.2), stringr (1.5.1), ggforce (0.4.2), grid (4.4.2), tmap (4.0), cowplot (1.1.3), scatterpie (0.2.5), gridExtra (2.3), ggrepel (0.9.6)

Python scripts were run in python version 3.12.7 and using the following packages: numpy (1.26.4), scikit-allel (1.3.13), pysam (0.23.3), plotly (5.24.1), scipy (1.14.1), matplotlib (3.9.3), pandas (2.2.3), tqdm (4.66.5)

## Bioinformatics analysis

Scripts can be run in the following order:
1. [batch_pipeline.sh](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/batch_pipeline.sh): mapping and variant calling of fastq files
2. [collect_vcf.sh](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/collect_vcf.sh): collect VCF files 
3. [create_merged_files.sh](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/create_merged_files.sh): create merged vcf
4. [realign-crt-haplotype.py](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/realign-crt-haplotype.py): re-align CRT haplotype
5. [analyse_vcf_EMAGEN.py](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/analyse_vcf_EMAGEN.py), [analyse_vcf_project.py](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/analyse_vcf_project.py): principal component analysis
6. [Haplotype_network_script.R](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/Haplotype_network_script.R): create haplotype networks
   
Prerequisite files
- bed file with genomic locations of amplicons can be found in input folder: [intervals.bed](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/input/intervals.bed)<br />
- Plasmodium falciparum reference file: [Pfalciparum.genome.fasta](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/input/Pfalciparum.genome.fasta)<br />
- Plasmodium falciparum gff file: [Pfalciparum.genome.modified.new.gff3](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/input/Pfalciparum.genome.modified.new.gff3)

## Statistical analysis

[Haplotype_corr_test.R](https://github.com/leenvh/EMAGEN/blob/main/statistical%20analysis/Haplotype_corr_test.R): haplotype correlation analysis<br />
[coocurance_analysis.R](https://github.com/leenvh/EMAGEN/blob/main/statistical%20analysis/coocurance_analysis.R): co-occurence analysis


## Visualisations

Graphs were made using R scripts and Jupyter notebook files.

- Figure 1: [Figure_1.R](https://github.com/leenvh/EMAGEN/blob/main/visualisations/Figure_1.R)
- Figure 2 and Figure 3A-B: [Figure_2_3A-B.R](https://github.com/leenvh/EMAGEN/blob/main/visualisations/Figure_2_3A-B.R)
- Figure 3C: [Upset_intersection_Plot.ipynb](https://github.com/leenvh/EMAGEN/blob/main/visualisations/Upset_intersection_Plot.ipynb)
- Figure 3D: [Co_occurrence_count_analysis.ipynb](https://github.com/leenvh/EMAGEN/blob/main/visualisations/Co_occurrence_count_analysis.ipynb)
- Figure 4 graphs are direct outputs from scripts in bioinformatics analysis folder
- Supplementary or Extended Data figures:
  * Haplotype correlation heatmap: [correlation_plot_heatmap.R](https://github.com/leenvh/EMAGEN/blob/main/visualisations/correlation_plot_heatmap.R)<br />
  * Co-occurence heatmap: [Co_occurence_plot_heatmap.R](https://github.com/leenvh/EMAGEN/blob/main/visualisations/Co_occurence_plot_heatmap.R)<br />
  * Heatmap NGS results: [MakingHeatmapOutline.ipynb](https://github.com/leenvh/EMAGEN/blob/main/visualisations/MakingHeatmapOutline.ipynb) <br />
  * Heatmap coverage data: [MakingCoverageHeatmapOutline.ipynb](https://github.com/leenvh/EMAGEN/blob/main/visualisations/MakingCoverageHeatmapOutline.ipynb) <br /> 
