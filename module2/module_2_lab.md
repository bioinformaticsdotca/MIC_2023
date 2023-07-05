---
layout: tutorial_page
permalink: /MIC_2023_Module2_lab
title: CBW-IMPACTT Microbiome Analysis 2023
header1: Workshop Pages for Students
header2: CBW-IMPACTT Microbiome Analysis Module 2 Lab
image: /site_images/CBW_Metagenome_icon.jpg
home: https://bioinformaticsdotca.github.io/
---

This tutorial is for the [2023 CBW-IMPACTT bioinformatics workshop](https://bioinformatics.ca/workshops-all/2023-cbw-impactt-microbiome-analysis/) running from July 5-7th with access to AWS servers. It is based on the Amplicon SOP v2 available on the [Microbiome Helper repository](https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2022.11)) and a [previous workshop designed by Monica Alvaro Fuss](https://github.com/LangilleLab/microbiome_helper/wiki/Microbiome-analysis-using-QIIME2-with-16S-data). 

**Authors:** Robyn Wright and Morgan Langille

## Table of Contents

1. First steps

2. Denoising the reads into amplicon sequence variants

3. Assign taxonomy to ASVs

4. Filtering resultant table

## Introduction

This lab provides a walkthrough of an end-to-end pipeline using the command line interface for the analysis of high-throughput marker gene data. Commonly used marker genes for microbiome analysis include the 16S ribosomal RNA (rRNA) for prokaryotes, 18S rRNA for eukaryotes, and the internal transcribed spacer (ITS) for fungi. In this tutorial, we will explore a 16S rRNA dataset from wild blueberry (*Vaccinium angustifolium*) soil communities from both bulk and rhizosphere soils associated with either natural or managed habitats. You can read more about the study in these papers:

- [Variation in Bacterial and Eukaryotic Communities Associated with Natural and Managed Wild Blueberry Habitats](https://apsjournals.apsnet.org/doi/10.1094/PBIOMES-03-17-0012-R)
- [Metagenomic Functional Shifts to Plant Induced Environmental Changes](https://www.frontiersin.org/articles/10.3389/fmicb.2019.01682/full#B50)

In this module (module 2) we will cover the basics of marker gene analysis from raw reads to filtered feature table. The pipeline described is embedded in the latest version of QIIME2 (Quantitative Insights into Microbial Ecology version 2023.2), which is a popular microbiome bioinformatics platform for microbial ecology built on user-made software packages called plugins that work on QIIME2 artifact or QZA files. Documentation for these plugins can be found in the [QIIME 2 user documentation](https://docs.qiime2.org/2023.2/), along with tutorials and other useful information. QIIME2 also provides interpretable visualizations that can be accessed by opening any generated QZV files within [QIIME2 View](https://view.qiime2.org/).
There are questions throughout this workshop that are aimed at ensuring you've understood the material - you can find the answers [here](https://github.com/bioinformaticsdotca/MIC_2023/blob/main/module2/module_2_answers.md) and are welcome to discuss the answers with us, but we will *not* be marking them.

#### Practical considerations

- If you run into problems copying the commands directly, try copying each line individually or typing them manually into the terminal. 

- Visualizing your data is always a good idea at any time during an analysis. In the interest of time, you will only be required to do so when you need to retrieve specific information from the visualization, but we have added some optional steps for you to explore if you wish.

- Due to the time and memory constraints of this workshop, we will provide some of the necessary output files.

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

## Getting set up

Open up your Terminal and log into the AWS server as you did previously.

QIIME2 is recommended to be run in an independent `conda` environment to ensure consistency between different versions of the packages it uses. For this instance, a QIIME2 environment was already installed, it just needs to be activated. [Here](https://docs.qiime2.org/2022.8/install/native/#install-qiime-2-within-a-conda-environment) are the instructions for installation if you need them. 

You can read more about conda/Anaconda [here](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html), but basically we often use it so that we can use *e.g.* program A requires version 2 of a program Z, that it is dependent on, while program B that we also have installed might require version 1 of program Z and installing version 2 might break program B. Conda makes sure that these separate versions stay separate from each another. It also often makes the installation of new packages very straightforward.

```
conda activate qiime2-2023.2
```

When you are finished this tutorial you can deactivate the conda environment using:

```
conda deactivate
```
*Note: don't run this yet! If you already did, re-run the above activate command.*

Create your tutorial directory, and then, move into this directory:

```
cd workspace
mkdir 16S_workshop/
```
```
cd 16S_workshop/
```

Link the tutorial data:

```
ln -s ~/CourseData/MIC_data/16S_data/raw_data/ .
ln -s ~/CourseData/MIC_data/16S_data/Blueberry_metadata_reduced.tsv .
```

If you log out of the server, make sure you change back into your directory and re-activate the conda environment before you resume the tutorial:

```
cd workspace/16S_workshop/
conda activate qiime2-2023.2
```

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

## 1. First steps

### 1.1. Inspect raw data

First, let's take a look at the directory containing our raw reads as well as our metadata file. 

```
ls raw_data
```
```
head Blueberry_metadata_reduced.tsv
```

**Question 1: How many samples are there?**

**Question 2: Into what group(s) are the samples classified?**

### *1.2. (Optional) Inspect raw read quality*

This can be important as a sanity check and to identify outlier samples with particularly low quality, read sizes, read depth, and other metrics. When dealing with raw reads the easiest method to use is a combination of [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/), as outlined in the [Microbiome Helper repository](https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2022.11)). Note that these tools are not packaged with QIIME 2 so you will need to install them separately.

### 1.3. Import FASTQs as QIIME2 artifact

To standardize QIIME 2 analyses and to keep track of provenance (i.e. a list of what commands were previously run to produce a file) a special format is used for all QIIME 2 input and output files called an "artifact" (with the extension QZA). The first step is to import the raw reads as a QZA file. We will first create a new directory.

```
mkdir reads_qza
```

```
qiime tools import \
    --type SampleData[PairedEndSequencesWithQuality] \
    --input-path raw_data/ \
    --output-path reads_qza/reads.qza \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt
```

All of the FASTQs are now in the single artifact file reads_qza/reads.qza. This file format can be a little confusing at first, but it is actually just a zipped folder. You can manipulate and explore these files better with the `qiime tools utilities` (e.g. `peek` and `view`).

### 1.4. Trim primers with cutadapt

Screen out reads that do not begin with primer sequence and remove primer sequence from reads using the [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html) QIIME 2 plugin. The below primers correspond to the 16S V6-V8 region (Bacteria-specific primer set).

```
qiime cutadapt trim-paired \
    --i-demultiplexed-sequences reads_qza/reads.qza \
    --p-cores 4 \
    --p-front-f ACGCGHNRAACCTTACC \
    --p-front-r ACGGGCRGTGWGTRCAA \
    --p-discard-untrimmed \
    --p-no-indels \
    --o-trimmed-sequences reads_qza/reads_trimmed.qza
```

**Question 3: What would happen if you ran this exact command on V4/V5-amplified sequences?**

### *1.4. (Optional) Summarize trimmed FASTQs*

Visualizing your output data is a good idea after any step. We still have a few preprocessing requirements to check off our list before denoising, so we can wait until these steps are complete to visualize our data. However, if you would like to see what paired-end reads look like before joining, run the following command and open the QZV file in [QIIME2 View](https://view.qiime2.org/):

```
qiime demux summarize \
    --i-data reads_qza/reads_trimmed.qza \
    --o-visualization reads_qza/reads_trimmed_summary.qzv
```

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

## 2. Denoising the reads into amplicon sequence variants

Different denoising tools require different levels of preprocessing before the actual denoising happens. For example, DADA2 performs read joining and quality filtering as part of the denoising step itself, and can be run directly after trimming the primers. Due to speed considerations, we'll be using Deblur instead, which requires that these steps be carried out separately. Guidelines for running DADA2 can be found [here](https://github.com/LangilleLab/microbiome_helper/wiki/QIIME2-DADA2-Quick-Reference). 

### 2.1. Join paired-end reads

Forward and reverse reads can be joined with [VSEARCH](https://github.com/torognes/vsearch) as shown below.

```
qiime vsearch merge-pairs \
    --i-demultiplexed-seqs reads_qza/reads_trimmed.qza \
    --o-merged-sequences reads_qza/reads_trimmed_joined.qza
```

### 2.2. Filter out low-quality reads

This command will filter out low-quality reads based on the default options.

```
qiime quality-filter q-score \
    --i-demux reads_qza/reads_trimmed_joined.qza \
    --o-filter-stats filt_stats.qza \
    --o-filtered-sequences reads_qza/reads_trimmed_joined_filt.qza
```

### 2.3. Summarize joined and filtered reads

It is a good idea at this point just to verify that there haven't been any substantial losses of reads, before going through the whole ASV process, at either the joining or quality-filtering steps above. You will also need to select a length to trim back to that maintains the largest/acceptable quantity of reads during denoising. 

```
qiime demux summarize \
    --i-data reads_qza/reads_trimmed_joined_filt.qza \
    --o-visualization reads_qza/reads_trimmed_joined_filt_summary.qzv
```

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

Now open the file in [QIIME2 View](https://view.qiime2.org/) and look at the Overview and Interactive Quality Plot tabs to explore your data and answer the following questions. 

**Question 4: How long are our forward reads? Why are there no reverse reads in our file?**

**Question 5: What would be a good trim length for our reads?**

### 2.4. Running Deblur

Running the [deblur](https://github.com/biocore/deblur) workflow will correct the raw reads into amplicon sequence variants (ASVs). This denoising tool filters out reads that either do match to known noise or that do not match with low similarity to the expected amplicon region. Note that the below command will retain singletons, which would have been filtered out unless we set --p-min-reads 1, and is for 16S sequences only. For other amplicon regions, you can either use the denoise-other option in the command and specify a reference database of sequences to use for positive filtering (as in the below versions for 18S and ITS) or use DADA2. We must also specify the trim length (previous step) to prevent getting an error message. 

```
qiime deblur denoise-16S \
    --i-demultiplexed-seqs reads_qza/reads_trimmed_joined_filt.qza \
    --p-trim-length 392 \
    --p-sample-stats \
    --p-jobs-to-start 4 \
    --p-min-reads 1 \
    --output-dir deblur_output
```
*Note: this command will probably take 10 minutes or so to run.*

When running 18S and ITS data with deblur you will need to specify custom reference files, which can be downloaded [here](http://kronos.pharmacology.dal.ca/public_files/MH/deblur_non16S_ref/). Further instructions on denoising 18S and ITS can be found on the [Microbiome Helper repository](https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2022.11))

### 2.5. Summarizing deblur output

Once a denoising pipeline has been run you can summarize the output table with the below command, which will create a visualization artifact for you to view. We will use this visualization later to determine the the cut-offs for filtering the table below, but for now you should mainly take a look at the visualization to ensure that sufficient reads have been retained after running deblur. This denoising tool filters out reads that either do match to known noise or that do not match with low similarity to the expected amplicon region. If your samples have very low depth after running deblur (compared to the input read depth) this could be a red flag that either you ran the tool incorrectly, you have a lot of noise in your data, or that deblur is inappropriate for your dataset.

```
qiime feature-table summarize \
    --i-table deblur_output/table.qza \
    --o-visualization deblur_output/deblur_table_summary.qzv
```

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

**Question 6: What is the mean sequencing depth per sample after denoising?**

**Question 7: Which sample has the least reads?**

## 3. Assign taxonomy to ASVs

You can assign taxonomy to your ASVs using a Naive-Bayes approach implemented in the [scikit learn](http://scikit-learn.org/stable/) Python library and the [SILVA](https://www.arb-silva.de/) or [UNITE](https://unite.ut.ee/) databases.

### *3.1. (Optional) Build or acquire taxonomic classifier*

This approach requires that a classifier be trained in advance on a reference database. We recommend users use a widely used classifier to help ensure there are no unexpected issues with the Naive-Bayes model. We previously maintained primer-specific classifiers, which theoretically can provide more accurate classifications, but we no longer do this due to concerns regarding issues with the trained models that are difficult to catch if only a couple people are running them.

The full-length 16S/18S classifier can be downloaded from the [QIIME 2 website](https://docs.qiime2.org/2022.11/data-resources/) (silva-138-99-nb-classifier.qza for the latest classifier).

Custom classifiers for the ITS region that we have generated from the UNITE database are available as well ([see downloads](http://kronos.pharmacology.dal.ca/public_files/MH/taxa_classifiers/qiime2-2020.8_classifiers) and [commands used to create these files](https://github.com/LangilleLab/microbiome_helper/wiki/Creating-QIIME-2-Taxonomic-Classifiers)):

-   Full ITS - fungi only (classifier_sh_refs_qiime_ver9_99_s_27.10.2022_ITS.qza)
-   Full ITS - all eukaryotes (classifier_sh_refs_qiime_ver9_99_s_all_27.10.2022_ITS.qza)

### 3.2. Run taxonomic classification

You can run the taxonomic classification with this command, which is one of the longest running and most memory-intensive command of the SOP. If you receive an error related to insufficient memory (and if you cannot increase your memory usage) then you can look into the --p-reads-per-batch option and set this to be lower than the default (which is dynamic depending on sample depth and the number of threads) and also try running the command with fewer jobs (e.g. set --p-n-jobs 1). **Because of the memory constraints of this workshop, you can copy the output into your folder:**

```
cp -r ~/CourseData/MIC_data/16S_data/taxa .
```

This is the command that you would usually run:
```
qiime feature-classifier classify-sklearn \
    --i-reads deblur_output/representative_sequences.qza \
    --i-classifier silva-138-99-nb-classifier.qza \
    --p-n-jobs 4 \
    --output-dir taxa
```
*Note: don't run this! If you started running it accidentally, you can stop it by pressing ```ctrl```+```c``` (at the same time). You can always use this to stop the command that is currently running.

### *3.3 (Optional) Assess subset of taxonomic assignments with BLAST*

As with all QZA files, you can export the output file to take a look at the classifications and confidence scores:

```
qiime tools export \
    --input-path taxa/classification.qza \
    --output-path taxa
```

The performance of the taxonomic classification is difficult to assess without a gold-standard reference, but nonetheless one basic sanity check is to compare the taxonomic assignments with the top BLASTn hits for certain ASVs. First, generate a QZV file for the denoised representative sequences in QIIME 2 by running:

```
qiime feature-table tabulate-seqs \
    --i-data deblur_output/representative_sequences.qza \
    --o-visualization deblur_output/representative_sequences.qzv
```

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

This QZV file tabulates the denoised sequences. Clicking on the nucleotide sequence links to a BLASTn search for that sequence.  By comparing these BLAST hits with the taxonomic assignment of ASVs generated above you can reassure yourself that the taxonomic assignments overall worked correctly. It's a good idea to select ~5 ASVs to BLAST for this validation, which should be from taxonomically different groups, such as different phyla, according to the taxonomic classifier.

## 4. Filtering resultant table

Filtering the denoised table is an important step of microbiome data analysis. You can see more details on this process in the [QIIME 2 filtering tutorial](https://docs.qiime2.org/2022.11/tutorials/filtering/).

### 4.1. Filter out rare ASVs

Based on the summary visualization created in step 2.5 above you can choose a cut-off for how frequent a variant needs to be (and optionally how many samples need to have the variant) for it to be retained. Here we will remove all ASVs that have a frequency of less than 0.1% of the mean sample depth. This cut-off excludes ASVs that are likely due to MiSeq bleed-through between runs (reported by Illumina to be 0.1% of reads). To calculate this cut-off you would identify the mean sample depth in the visualization created in step 2.5, multiply it by 0.001, and round to the nearest integer.

Here I've added in the number that you should use (11), so that these files are the same for everyone for some of the next steps, but make sure that you're happy with how this is determined:
```
qiime feature-table filter-features \
    --i-table deblur_output/table.qza \
    --p-min-frequency 11 \
    --p-min-samples 1 \
    --o-filtered-table deblur_output/deblur_table_filt.qza
```

### 4.2. Filter out contaminant and unclassified ASVs

Once we have assigned taxonomy to our ASVs we can use that information to remove ASVs which are likely contaminants or noise based on the taxonomic labels. Two common contaminants in 16S sequencing data are mitochondrial and chloroplast 16S sequences, which can be removed by excluding any ASV which contains those terms in its taxonomic label. It can also be **sometimes** useful to exclude any ASV that is unclassified at the phylum level since these sequences could be noise (e.g. possible chimeric sequences). Note that if your data has not been classified against the default database you may need to change p__ to be a string that enables phylum-level assignments to be identified or simply omit that line.

In general though, it can be very informative if your sequencing reads are coming back with significant amounts of unclassified ASVs as it can indicate upstream analysis problems or indicate you are studying a poorly characterized environment where you have a good chance of identifying a lot of novel phyla. Therefore, our recommendation is to not filter out the unclassified sequences by default:

If you do want to exclude all the unclassifieds, then run the following instead (with the above caveat that the line might need to be adapted according to the tax database used):

```
qiime taxa filter-table \
    --i-table deblur_output/deblur_table_filt.qza \
    --i-taxonomy taxa/classification.qza \
    --p-include p__ --p-exclude mitochondria,chloroplast \
    --o-filtered-table deblur_output/deblur_table_filt_contam.qza
```

### *4.3. (Optional) Exclude low-depth samples*

Often certain samples will have quite low depth after these filtering steps, which can be excluded from downstream analyses since they will largely add noise. There is no single cut-off that works best for all datasets, but researchers often use minimum cut-offs within the range of 1000 to 4000 reads. You can also use a cut-off much lower than this if you want to retain all samples except those that failed entirely (e.g. depth < 50 reads). 

Ideally you would choose this cut-off after visualizing rarefaction curves to determine at what read depth the richness of your samples plateaus and choose a cut-off as close to this plateau as possible while retaining sufficient sample size for your analyses. We learn more about rarefaction curves in part two.

### 4.4. Subset and summarize filtered table

Check output after filtering. 

```
qiime feature-table summarize \
    --i-table deblur_output/deblur_table_filt_contam.qza \
    --o-visualization deblur_output/deblur_table_filt_contam_summary.qzv
```

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**

**Question 8: What is the maximum sequencing depth across all samples?**

Happy? Copy a final table. 

```
cp deblur_output/deblur_table_filt_contam.qza deblur_output/deblur_table_final.qza
```

Once we have our final filtered table we will need to subset the QZA of ASV sequences to the same set. You can exclude any removed ASVs from the sequence file with this command:

```
qiime feature-table filter-seqs \
    --i-data deblur_output/representative_sequences.qza \
    --i-table deblur_output/deblur_table_final.qza  \
    --o-filtered-data deblur_output/rep_seqs_final.qza
```

Finally, you can make a new summary of the final filtered abundance table:

```
qiime feature-table summarize \
    --i-table deblur_output/deblur_table_final.qza \
    --o-visualization deblur_output/deblur_table_final_summary.qzv
```

We will also usually want to export this final table as well as a fasta file of the final ASVs, which we can do with the following commands:
```
#first export the taxonomic classifications as a table
qiime tools export \
  --input-path taxa/classification.qza \
  --output-path taxa

#rename some of this table
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxa/taxonomy.tsv

#export the final table (by default this will be in .biom format) - we call this a feature table, and within it, 
#rows will be ASVs and columns will be samples, with each cell containing the number of reads assigned to an ASV within that sample
qiime tools export \
  --input-path deblur_output/deblur_table_final.qza \
  --output-path exports

#add the taxonomy information to the feature table
biom add-metadata \
  -i exports/feature-table.biom \
  -o exports/feature-table_w_tax.biom \
  --observation-metadata-fp taxa/taxonomy.tsv \
  --sc-separated taxonomy
 
#convert the .biom feature table with taxonomy into a .txt file
biom convert \
  -i exports/feature-table_w_tax.biom \
  -o exports/feature-table_w_tax.txt \
  --to-tsv \
  --header-key taxonomy

#export the ASV sequences as a fasta file
qiime tools export \
  --input-path deblur_output/rep_seqs_final.qza \
  --output-path exports
```
You can go ahead an open this file up in Excel if you like. We'll be doing some visualisations within QIIME2 in the next workshop, but you could also analyse these files in R, Python, etc., or in web-based platforms like [Microbiome Analyst](https://www.microbiomeanalyst.ca/).

**Remember that you can look at all of the files in your ```workspace``` folder by going to http://##.uhn-hpc.ca/ (you might need to refresh the page regularly!)**
