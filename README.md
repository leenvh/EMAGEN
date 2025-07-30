
This repository contains code used to process, analyze and visualise Plasmodium falciparum genomic data from raw fastq files as described in the manuscript [Convergent Evolution of Artemisinin and Chloroquine Resistance in Ethiopian Plasmodium falciparum Parasites](https://verixiv.org/articles/2-162)

## Bioinformatics analysis

Bash scripts scripts rely on functioning installations of bwa, python, samtools, gatk4, bcftools, and snpEff.

Scripts can be run in the following order:
1. [batch_pipeline.sh](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/batch_pipeline.sh): mapping and variant calling of fastq files
2. [collect_vcf.sh](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/collect_vcf.sh): collect VCF files 
3. [create_merged_files.sh](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/create_merged_files.sh): create merged vcf
4. [realign-crt-haplotype.py](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/realign-crt-haplotype.py): re-align CRT haplotype
5. [analyse_vcf_EMAGEN.py](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/analyse_vcf_EMAGEN.py), [analyse_vcf_project.py](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/analyse_vcf_project.py): principal component analysis
6. [Haplotype_network_script.R](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/Haplotype_network_script.R): create haplotype networks
   
Prerequisite files
- bed file with genomic locations of amplicons can be found in input folder: [intervals.bed](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/input/intervals.bed)<br />
- Plasmodium falciparum reference file: [Pfalciparum.genome.fasta](https://github.com/leenvh/EMAGEN/blob/main/bioinformatics%20analysis/input/Pfalciparum.genome.fasta)

## Statistical analysis

[Haplotype_corr_test.R](https://github.com/leenvh/EMAGEN/blob/main/Haplotype_corr_test.R): haplotype correlation analysis


## Visualisations

Graphs were made using R scripts and Jupyter notebook files.

- Figure 1: [Figure_1.R](https://github.com/leenvh/EMAGEN/blob/main/visualisations/Figure_1.R)
- Figure 2 and Figure 3A-B: [Figure_2_3A-B.R](https://github.com/leenvh/EMAGEN/blob/main/visualisations/Figure_2_3A-B.R)
- Figure 3C: [Upset_intersection_Plot.ipynb](https://github.com/leenvh/EMAGEN/blob/main/visualisations/Upset_intersection_Plot.ipynb)
- Figure 3D: [Co_occurrence_count_analysis.ipynb](https://github.com/leenvh/EMAGEN/blob/main/visualisations/Co_occurrence_count_analysis.ipynb)
- Figure 4 graphs are direct outputs from scripts in bioinformatics analysis folder
- Supp figures:
  * Haplotype correlation heatmap: [correlation_plot_heatmap.R](https://github.com/leenvh/EMAGEN/blob/main/visualisations/correlation_plot_heatmap.R)
  * Heatmap NGS results: [MakingHeatmapOutline.ipynb](https://github.com/leenvh/EMAGEN/blob/main/visualisations/MakingHeatmapOutline.ipynb) <br />
  * Heatmap coverage data: [MakingCoverageHeatmapOutline.ipynb](https://github.com/leenvh/EMAGEN/blob/main/visualisations/MakingCoverageHeatmapOutline.ipynb) <br /> 
