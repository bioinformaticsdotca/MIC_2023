---
layout: tutorial_page
permalink: /MIC_2023_Module2_answers
title: CBW-IMPACTT Microbiome Analysis 2023
header1: Workshop Pages for Students
header2: CBW-IMPACTT Microbiome Analysis Module 2 Lab Answers
image: /site_images/CBW_Metagenome_icon.jpg
home: https://bioinformaticsdotca.github.io/
---


**Question 1: How many samples are there?**

We have 12 samples. The raw data folder contains one fastq.gz file for each set of the forward and reverse sequences (labelled R1 and R2, respectively) for each of the samples (BBxxx).

**Question 2: Into what group(s) are the samples classified?**

The samples are classified into two different groups: either bulk or rhizosphere soil and forest or managed sites.

**Question 3: What would happen if you ran this exact command on V4/V5-amplified sequences?**

```Cutadapt``` will only trim reads that match the specified primer sequence. Therefore, most reads would be discarded because we are including the ```--p-discard-untrimmed``` option.

**Optional Question: How long are the reads now? Is this what you were expecting?**

You should see that the reads are around 280 bp long now, because the ~20 bp primers have been trimmed from the ends and they were ~300 bp to start). 

**Question 4: How long are our forward reads? Why are there no reverse reads in our file?**

The forward read median length is 405 nucleotides. There are no reverse reads because forward and reverse reads were merged into one sequence during read joining.

**Question 5: What would be a good trim length for our reads?**

There is no one right answer for this question, but a trim length of 392 nucleotides is good enough to maintain a proper balance between depth and quality.

**Question 6: What is the mean sequencing depth per sample after denoising?**

The mean sequencing depth (frequency) across all denoised samples is 11,065.25 reads.

**Question 7: Which sample has the least reads?**

Sample BB203 has the lowest sequencing depth (7,593 reads).

**Question 8: What is the maximum sequencing depth across all samples?**

The final maximum sequencing depth is 11,690 reads.
