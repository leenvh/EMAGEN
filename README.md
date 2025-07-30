
This repository contains code used to process, analyze and visualise Plasmodium falciparum genomic data from raw fastq files as described in the manuscript [Convergent Evolution of Artemisinin and Chloroquine Resistance in Ethiopian Plasmodium falciparum Parasites](https://verixiv.org/articles/2-162)

## Bioinformatics analysis


Bed file with genomic locations of amplicons: [intervals.bed](https://github.com/leenvh/EMAGEN/blob/main/intervals.bed)<br />
Mapping and variant calling of fastq files: [batch_pipeline.sh](https://github.com/leenvh/EMAGEN/blob/main/batch_pipeline.sh), [run_variant_calling.sh](https://github.com/leenvh/EMAGEN/blob/main/run_variant_calling.sh)<br />
Create merged vcf: [collect_vcf.sh](https://github.com/leenvh/EMAGEN/blob/main/collect_vcf.sh), [create_merged_files.sh](https://github.com/leenvh/EMAGEN/blob/main/create_merged_files.sh)<br />
Re-align CRT haplotype: [realign-crt-haplotype.py](https://github.com/leenvh/EMAGEN/blob/main/realign-crt-haplotype.py)<br />
Principal component analysis: [analyse_vcf_EMAGEN.py](https://github.com/leenvh/EMAGEN/blob/main/analyse_vcf_EMAGEN.py), [analyse_vcf_project.py](https://github.com/leenvh/EMAGEN/blob/main/analyse_vcf_project.py)<br />
Haplotype network: [Haplotype_network_script.R](https://github.com/leenvh/EMAGEN/blob/main/Haplotype_network_script.R)<br />
Haplotype correlation analysis: [Haplotype_corr_test.R](https://github.com/leenvh/EMAGEN/blob/main/Haplotype_corr_test.R) <br />
Haplotype correlation heatmap: [correlation_plot_heatmap.R](https://github.com/leenvh/EMAGEN/blob/main/correlation_plot_heatmap.R) with [corr_data_output.rds](https://github.com/leenvh/EMAGEN/blob/main/corr_data_output.rds) <br /> 
Heatmap NGS results: [MakingHeatmapOutline.ipynb](https://github.com/leenvh/EMAGEN/blob/main/MakingHeatmapOutline.ipynb) <br /> 
Heatmap coverage data: [MakingCoverageHeatmapOutline.ipynb](https://github.com/leenvh/EMAGEN/blob/main/MakingCoverageHeatmapOutline.ipynb) <br /> 
