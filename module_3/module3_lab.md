---
layout: tutorial_page
permalink: /MIC_2023_Module3_lab
title: CBW-IMPACTT Microbiome Analysis 2023
header1: Workshop Pages for Students
header2: CBW-IMPACTT Microbiome Analysis Module 3 Lab
image: /site_images/CBW_Metagenome_icon.jpg
home: https://bioinformaticsdotca.github.io/
--- 

Microbial diversity metrics and data visualization with QIIME2


This tutorial is for the 2023 CBW IMPACTT workshop running from July 5th-7th with access to AWS instance. It is based on the Amplicon SOP v2 available on the [Microbiome Helper repository](https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2022.11)) and previous workshops designed by Diana Haider and Robert Beiko and Monica Alvaro Fuss and Morgan Langille

Authors: Robyn Wright & Hena Ramay

# Table of Contents
1. Build tree
2. Generate rarefaction curves
3. Calculating diversity metrics and generating ordination plots
4. Generate stacked bar chart of taxa relative abundances
5. Identifying differentially abundant features with ANCOM
6. Exporting the final abundance, profile and sequence files

# Introduction

This lab provides a walkthrough of an end-to-end pipeline using the command line interface for the analysis of high-throughput marker gene data. Commonly used marker genes for microbiome analysis include the 16S ribosomal RNA (rRNA) for prokaryotes, 18S rRNA for eukaryotes, and the internal transcribed spacer (ITS) for fungi. In this tutorial, we will explore a 16S rRNA dataset from wild blueberry (Vaccinium angustifolium) soil communities from both bulk and rhizosphere soils associated with either natural or managed habitats. You can read more about the study in these papers:

* [Variation in Bacterial and Eukaryotic Communities Associated with Natural and Managed Wild Blueberry Habitats](https://apsjournals.apsnet.org/doi/10.1094/PBIOMES-03-17-0012-R)

* [Metagenomic Functional Shifts to Plant Induced Environmental Changes](https://www.frontiersin.org/articles/10.3389/fmicb.2019.01682/full#B50)



In part 1 we covered the basics of marker gene analysis from raw reads to filtered feature table, and in part 2 we will explore examples of downstream analyses that can be used to draw biological conclusions from these data. The pipeline described is embedded in the latest version of QIIME2 (Quantitative Insights into Microbial Ecology version 2023.2), which is a popular microbiome bioinformatics platform for microbial ecology built on user-made software packages called plugins that work on QIIME2 artifact or QZA files. Documentation for these plugins can be found in the [QIIME 2 user documentation](https://docs.qiime2.org/2023.2/), along with tutorials and other useful information. QIIME2 also provides interpretable visualizations that can be accessed by opening any generated QZV files within [QIIME2 View](https://view.qiime2.org/).

There are questions throughout this workshop that are aimed at ensuring you've understood the material - you can find the answers [here]() and are welcome to discuss the answers with us, but we will not be marking them.

### Practical considerations

If you run into problems copying the commands directly, try copying each line individually or typing them manually into the terminal.

Visualizing your data is always a good idea at any time during an analysis. In the interest of time, you will only be required to do so when you need to retrieve specific information from the visualization, but we have added some optional steps for you to explore if you wish.

If you run out of memory, we will provide the necessary output files.

### Activate QIIME2 conda environment

QIIME2 is recommended to be run in an independent conda environment to ensure consistency between different versions of the packages it uses. For this instance, a QIIME2 environment was already installed, it just needs to be activated. [Here](https://docs.qiime2.org/2022.8/install/native/#install-qiime-2-within-a-conda-environment) are the instructions for installation if you need them.

```
conda activate qiime2-2023.2
```

When you are finished this tutorial you can deactivate the conda environment using:

```
conda deactivate
```

## 1. Build tree with [SEPP QIIME 2 plugin](https://github.com/qiime2/q2-fragment-insertion)
[SEPP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5904434/) is one method for placing short sequences into a reference phylogenetic tree. This is a useful way of determining a phylogenetic tree for your ASVs. For 16S data you can do this with q2-fragment-insertion using the below command. Again, due to memory constraints, you can copy the output into your folder with the following command:

```
cp ~/CourseData/MIC_data/16S_data/asvs-tree.qza .
```
```
cp ~/CourseData/MIC_data/16S_data/insertion-placements.qza .
```

DONOT RUN THIS COMMAND
```
qiime fragment-insertion sepp \
    --i-representative-sequences deblur_output/rep_seqs_final.qza \
    --i-reference-database /home/shared/microbiome_amplicon/sepp-refs-gg-13-8.qza \
    --o-tree asvs-tree.qza \
    --o-placements insertion-placements.qza \
    --p-threads 4
```

Note that if you do not already have this file locally you will need to download sepp-refs-gg-13-8.qza as specified in the [fragment-insertion instructions](https://library.qiime2.org/plugins/q2-fragment-insertion/16/). You can specify custom reference files to place other amplicons, but the easiest approach for 18S and ITS data is to instead create a de novo tree as outlined in the [Microbiome Helper repository](https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2022.11)).

## 2. Generate rarefaction curves

A key quality control step is to plot rarefaction curves for all of your samples to determine if you performed sufficient sequencing. The below command will generate these plots (make sure you have the correct maximum sequencing depth as per your filtered feature table).

```
qiime diversity alpha-rarefaction \
    --i-table deblur_output/deblur_table_final.qza \
    --p-max-depth 11536 \
    --p-steps 20 \
    --i-phylogeny asvs-tree.qza \
    --o-visualization rarefaction_curves.qzv
```

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

**Question 1: What is a good rarefaction depth for diversity analysis?**



## 3. Calculating diversity metrics and generating ordination plots

Common **alpha** and **beta-diversity** metrics can be calculated with a single command in QIIME 2. In addition, ordination plots (such as PCoA plots for weighted UniFrac distances) will be generated automatically as well. This command will also rarefy all samples to the sample sequencing depth before calculating these metrics ( is a placeholder for the lowest reasonable sample depth; samples with depth below this cut-off will be excluded). In our case X can be 4000 to make sure that we have all the samples included. 

```
qiime diversity core-metrics-phylogenetic \
    --i-table deblur_output/deblur_table_final.qza \
    --i-phylogeny asvs-tree.qza \
    --p-sampling-depth X  \
    --m-metadata-file Blueberry_metadata_reduced.tsv \
    --p-n-jobs-or-threads 4 \
    --output-dir diversity
```

While `qiime diversity core-metrics-phylogenetic` forces you to specify sampling depth hence makes you rarefy samples, if you use a diversity 
measure directly you don't have to specify the depth. There is a lot of debate on rarefying or not, as general rule, you can skip it if the sample depth is within 10x of the smallest sample depth. 

### Alpha diversity visualization and significance test

For alpha diversity visualizations, you will need to produce boxplots comparing the different categories in your metadata file. For example, to create boxplots comparing the Shannon alpha-diversity metric you can use this command:

```
qiime diversity alpha-group-significance \
    --i-alpha-diversity diversity/shannon_vector.qza \
    --m-metadata-file Blueberry_metadata_reduced.tsv \
    --o-visualization diversity/shannon_compare_groups.qzv
```
**Hint: you will need to change this command for the other alpha diversity metrics. You can see the other metrics available by running ls diversity/*_vector.qza**

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

**Note that you can also export (see below) this or any other diversity metric file (ending in .qza) and analyze them with a different program.**

**Question 2: are there any significant differences in alpha diversity between any of our metadata categories?**

**Question 3: which metadata category appears to provide more separation in the beta diversity PCoA plots?**

### Beta diversity visualization and significance test

Here we are going to use PERMANOVA from the beta-group-significance function to see if out groups in the Description_3 variable are significantly different or not. By default the function performs PERMANOVA method.The PERMANOVA implementation here is one-way. To include more than one variable with potential interactions, use `qiime diversity adonis`

``` 
qiime diversity beta-group-significance \
    --i-distance-matrix diversity/bray_curtis_distance_matrix.qza \
    --m-metadata-file Blueberry_metadata_reduced.tsv \
    --m-metadata-column Description_3 \
    --o-visualization beta_group_sig_permanova
```

It is recommended to test for dispersion of groups especially if PERMANOVA shows signifcant difference in groups. This is to test if dispersions of the groups being compared is different or not. If the dispersion test has a significant p-value that means that group dispersions are different and you should be careful with reporting your PERMANOVA results.

```
qiime diversity beta-group-significance \
    --i-distance-matrix diversity/bray_curtis_distance_matrix.qza \
    --m-metadata-file Blueberry_metadata_reduced.tsv \
    --m-metadata-column Description_3 \
    --o-visualization beta_group_sig_permdisp \
    --p-method "permdisp"
```


## 4. Generate stacked bar chart of taxa relative abundances

Another useful output is the interactive stacked bar-charts of the taxonomic abundances across samples, which can be output with this command:
```
qiime taxa barplot \
    --i-table deblur_output/deblur_table_final.qza \
    --i-taxonomy taxa/classification.qza \
    --m-metadata-file Blueberry_metadata_reduced.tsv \
    --o-visualization taxa/taxa_barplot.qzv
```
**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

**Question 4: can you identify any patterns between the metadata groups?**

## 5. Identifying differentially abundant features with [ANCOM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450248/)

ANCOM is one method to test for differences in the relative abundance of features between sample groupings. It is a compositional approach that makes no assumptions about feature distributions. However, it requires that all features have non-zero abundances so a pseudocount first needs to be added (1 is a typical pseudocount choice):

```
qiime composition add-pseudocount \
    --i-table deblur_output/deblur_table_final.qza \
    --p-pseudocount 1 \
    --o-composition-table deblur_output/deblur_table_final_pseudocount.qza
```

Then ANCOM can be run with this command; note that CATEGORY is a placeholder for the text label of your category of interest from the metadata file:

```
qiime composition ancom \
    --i-table deblur_output/deblur_table_final_pseudocount.qza \
    --m-metadata-file Blueberry_metadata_reduced.tsv \
    --m-metadata-column CATEGORY \
    --output-dir ancom_output
```

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

**Question 5: Does ANCOM identify any differentially abundant taxa between any of the metadata groups? If so, which one(s)?**

## 6. Exporting the final abundance, profile and sequence files
If you didn't do this at the end of Module 2, to get the BIOM file (with associated taxonomy) and FASTA file (one per ASV) for your dataset to plug into other programs you can use the commands below.

To export the FASTA:

```
qiime tools export --input-path deblur_output/rep_seqs_final.qza --output-path deblur_output_exported
```

To export the BIOM table (with taxonomy added as metadata):
```
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxa/taxonomy.tsv
qiime tools export --input-path deblur_output/deblur_table_final.qza --output-path deblur_output_exported

biom add-metadata -i deblur_output_exported/feature-table.biom -o deblur_output_exported/feature-table_w_tax.biom --observation-metadata-fp taxa/taxonomy.tsv --sc-separated taxonomy

biom convert -i deblur_output_exported/feature-table_w_tax.biom -o deblur_output_exported/feature-table_w_tax.txt --to-tsv --header-key taxonomy
```
