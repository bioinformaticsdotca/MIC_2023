---
layout: tutorial_page
permalink: /MIC_2023_Module4_lab
title: CBW-IMPACTT Microbiome Analysis 2023
header1: Workshop Pages for Students
header2: CBW-IMPACTT Microbiome Analysis Module 4 Lab
image: /site_images/CBW_Metagenome_icon.jpg
home: https://bioinformaticsdotca.github.io/
---


This tutorial is for the [2023 CBW-IMPACTT bioinformatics workshop](https://bioinformatics.ca/workshops-all/2023-cbw-impactt-microbiome-analysis/) running in Calgary from July 5-7th. This tutorial has been adapted from the [2022 IMPACTT bioinformatics workshops](https://github.com/LangilleLab/microbiome_helper/wiki/IMPACTT-workshop-tutorials-2022), but it now uses different samples that match those in the [amplicon sequencing part of this workshop](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-IMPACTT-2023-Microbiome-analysis-using-QIIME2-with-16S-data). We have also added the part that includes visualisation of results in Phyloseq, and removed parts related to initial quality control of reads.

**Authors**: [Robyn Wright](mailto:robyn.wright@dal.ca) and Morgan Langille

### Table of Contents

#### Taxonomic annotation

1. Initial setup

2. Generating Taxonomic Profiles with Kraken2

3. Running Bracken

4. Visualising the results using Phyloseq in R

#### Functional annotation

1. Determining functional profiles using MMseqs2

2. Adding descriptions and stratifying the file

3. Visualising with JARRVIS

### Bioinformatic tool citations

Properly citing bioinformatic tools is important to ensure that developers are getting the proper recognition. If you use the commands in this tutorial for your own work be sure to cite the relevant tools listed below!

* **GNU Parallel** ([website](https://www.gnu.org/software/parallel/), [paper](https://www.usenix.org/system/files/login/articles/105438-Tange.pdf))
* **Microbiome Helper** ([website](https://github.com/LangilleLab/microbiome_helper/wiki), [paper](http://msystems.asm.org/content/2/1/e00127-16))
* **Kraken2** ([website](https://ccb.jhu.edu/software/kraken2/), [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0))
* **Bracken** ([website](https://ccb.jhu.edu/software/bracken/), [paper](https://peerj.com/articles/cs-104/))
* **Phyloseq** ([website](https://joey711.github.io/phyloseq/), [paper](https://doi.org/10.1371/journal.pone.0061217))
* **Vegan** ([website](https://vegandevs.github.io/vegan/))
* **MMseqs2** ([website](https://github.com/soedinglab/mmseqs2), [paper](https://www.nature.com/articles/nbt.3988))
* **JARRVIS** ([website](https://github.com/dhwanidesai/JarrVis))

## Taxonomic annotation introduction

The main goal of this tutorial is to introduce students to different approaches for the taxonomic profiling of metagenomic data from microbiome samples. We want to emphasize that there is not a one-size-fits-all pipeline for analyzing MGS data. For example some individuals may opt to examine their data using a marker-based approach such as [MetaPhlAn](https://huttenhower.sph.harvard.edu/metaphlan/) while others may opt to use a kmer based strategy such as Kraken2+Bracken. Furthermore, in a subsequent module in this workshop, you will also learn about another approach for examining microbiome data using metagenomic assembly and binning.

Throughout this tutorial there will be questions are aimed at helping students understand the various steps in reference based metagenomic profiling. [The answers are found on this page]().

This tutorial is based on our Microbiome Helper [Metagenomics SOP](https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Standard-Operating-Procedure-v3) and we are also working on a list of resources and papers that may be useful for beginners (or those that are not beginners!) in the microbiome field [here](https://github.com/LangilleLab/microbiome_helper/wiki/Microbiome-for-beginners).

## Initial setup

Usually the first steps in metagenome analysis would involve: 
* (1) Quality control - checking that all of the reads are of sufficient quality to use in our analysis - and primer trimming - removal of any of the primers/adapters used in sequencing
* (2) removal of unwanted reads - we usually check that there are no PhiX reads left over from the sequencing (usually added into Illumina sequencing runs for internal quality control) and may also remove any host-associated reads, for example, in human studies we would usually remove reads that map to the human genome
* (3) joining of the forward and reverse reads

Both of the first two steps can be performed using [Kneaddata](https://huttenhower.sph.harvard.edu/kneaddata/). The third step could be with a program that checks for matches between the reads, but as metagenome reads often don't have a large - or any - overlap, we may just simply concatenate the reads together. We've chosen to focus more on the taxonomic/functional annotations/downstream interpretation of our data in todays workshop, but you can read/work through the [previous workshop that we ran](https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-(taxonomic-annotation;-IMPACTT-2022)) for a complete overview of of these preliminary steps that you should definitely run if you have your own data.

First open up your terminal and log in to the server instance.

Next, we'll make a new directory to work from, change into it and then create a link to the sequence reads we'll be using. 
```
cd workspace/
mkdir metagenome_workshop
cd metagenome_workshop/
ln -s ~/CourseData/MIC_data/metagenome_data/cat_reads/ .
```
*Note: the ```ln -s``` command is an alternative to copying a file across and is like creating a desktop shortcut to a folder; we don't need to take up additional storage space but will still be able to use these files. It won't always be appropriate if we need to make modifications to files, but works well for cases like this where we wouldn't want to modify our raw data anyway.*

If you want to take a look at the read files, you can do this with *e.g.* ```less cat_reads/BB190.fastq``` (press ```q``` to get out of this view) of ```head cat_reads/BB190.fastq```

Next, we'll activate the conda environment that we'll be using for this part of the workshop:
```
conda activate taxonomic
```

If, during the workshop, you get logged out, you can get back to where you need to be by running:
```
cd workspace/metagenome_workshop/
conda activate taxonomic
```

## Generating taxonomic profiles with Kraken2

After we have performed quality control on all of our reads the next step is to use [Kraken2](https://ccb.jhu.edu/software/kraken2/) to classify them with taxonomic labels. Kraken2 is one of many popular tools for classifying metagenomic reads into taxonomic profiles. In this case Kraken2 uses a kmer approach to assign reads to specific taxonomic lineages (look back at the lecture slides to remind yourself how kmers are used to classify reads).

There are several choices that we can make when classifying our samples that can have a very significant impact on both the proportion of reads in samples that we are able to classify as well as the taxa that are classified (if a taxon doesn't exist in our database, then we won't be able to classify it in our samples). The biggest choices here are the reference database that we use as well as the confidence threshold. We have recently published a [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000949) that addresses this in some detail, but we'll summarise some of the most important points below.

As it will take about ten minutes for Kraken2 to run, we will just start running the commands and then there are a few things that we should understand about what we're doing.

First of all, we're going to make some directories to hold our results:

```
mkdir kraken2_kreport
mkdir kraken2_outraw
```

When you're running Kraken2 you should always be careful that you've created these folders before you run it! It won't give you an error message telling you that the folders don't exist, but it won't be able to create any output files so your time running it will have been wasted if you don't have the folders made. 

Next we will run the Kraken2 command using parallel:
```
parallel -j 1 --eta --dry-run 'kraken2 --db  ~/CourseData/MIC_data/metagenome_data/k2_pluspf_08gb_20230314/ --threads 4 --output kraken2_outraw/{1/.}_{2}_minikraken.kraken.txt --report kraken2_kreport/{1/.}_{2}_minikraken.kreport --use-names {1} --confidence {2}' ::: cat_reads/*.fastq ::: 0.0 0.1
```

Notice here that we now have two arguments going into parallel with ```:::``` at the end, so they are named ```{1}``` and ```{2}``` inside the command. Parallel is a really useful tool that allows us to run the same command on multiple files and with multiple parameters at once. In this workshop, we haven't given a tutorial on using parallel, but if you'd like to see that, it's in [this previous workshop](https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-(taxonomic-annotation;-IMPACTT-2022)#1-introducing-gnu-parallel). We also have some special notation with parallel where it will use {1} or {2} inside the command, and we can additionally use {1/} to remove the folder name from this file name, or {1.} to remove the extension (or {1/.} to remove both). In our example, the first file for {1} could be ```cat_reads/BB190.fastq```, and {1/} would be ```BB190.fastq```, {1.} would be ```cat_reads/BB190``` and {1/.} would be ```BB190```.

Check to make sure the command you are running is the one you expected using the --dry-run flag. Once you are sure remove the flag and run the command. This is really useful when you have a command that will take a really long time to run, so you don't waste time running the wrong command.

There are a lot of different options that we've specified in the kraken2 command:
* ```parallel -j 1 --eta --dry-run``` - here we've specified the number of jobs for parallel to use (```-j 1```), we've asked parallel to tell us a rough time to completion (```--eta```) and we're checking the command that will be run before running it (```--dry-run```).
* ```--db``` - This option points to the location of the database you want Kraken2 to use for classification. We have already set this AWS instance with the minikraken database.
* ```--threads 4```-This option indicates that we want to use four processors for each Kraken2 job
* ```--output``` - This argument points to the file name we want to output our classifications to.
* ```--report``` - This argument points to the file name we want to output a more detailed classifcation report for each sample. These more in-depth classification reports are required to correct abundance estimates using bracken. 
* ```--use-names``` - This indicates we want the classifier to output the taxon names rather than their NCBI IDs. 
* ```--confidence``` - This allows us to set the confidence threshold so as reads with a low proportion of their k-mers classified as the taxon will be unclassified in the output (see below).

Some summary information will be output to the console when you run this command including the number of sequences that were classified versus unclassified in each sample.

### Databases

Unlike metagenomic assembly, we will be classifying our sequenced reads based on how similar they are to reads we have already classified in a reference database. In the [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual), the Kraken2 developers recommend building a "standard" database that is comprised of the bacterial, archaeal and viral domains as well as the human genome and a collection of known vectors (UniVec_Core). There are commands given in the manual for doing this, and a relatively up-to-date version of this database is also available for download via the Langmead lab [here](https://benlangmead.github.io/aws-indexes/k2). This database requires approximately 100 GB RAM/memory to build, and 40 GB RAM to run. The Langmead Lab also provides some secondary databases that are capped at different sizes, and also that contain other domains in addition to the bacteria, archaea and viruses. In our paper, we found that the best database to use with a range of simulated and mock samples (where the real composition of the sample was known) was one that contained all genomes from the [NCBI Reference Sequence Database](https://www.ncbi.nlm.nih.gov/refseq/) (NCBI RefSeq), and can be downloaded from [here](https://www.dropbox.com/sh/lvlz2wpsssvsrad/AAC-BkJja8LvlDoNDB4qgnHNa?dl=0). However, this database is approximately 1.2 Tb and will therefore be too large for many researchers to use. We found that bigger databases weren't always better, and without access to a server with very large amounts of RAM, the best options were smaller, curated databases - we recommend that everyone thinks carefully about this choice prior to use on their own samples, and provided some guidance on this in the discussion of the paper.

Because we only have a limited amount of RAM/memory on our instances here, we will be using the standard database capped at 8 Gb, which we'll call the minikraken database. The smaller size means that we can run it on workstations with a lower amount of RAM, but this comes at the cost of reduced performance. Also note that by default, kraken uses NCBI taxonomy (and requires a taxonomy structure similar to that of the NCBI taxonomy, where all taxa at all ranks are assigned a taxonomic ID, and the "parent" taxon is known for all), however, there are a number of other databases that users have made using different styles of taxonomy such as the [Genome Taxonomy Database](https://gtdb.ecogenomic.org/).

### Confidence thresholds

As mentioned above, Kraken2, uses exact alignment of k-mers (sub-sequences of length k) to a reference database for query read classification. Because k-mers within a read could map to multiple sequences within the database when using Kraken2, a lowest common ancestor (LCA) approach is used to classify the sequence, and the number of k-mers mapping to a given taxon are calculated as a proportion of the number of k-mers within that sequence that do not have ambiguous nucleotides. This proportion of k-mers mapping to a taxon is called the “confidence”, and there is a confidence threshold that can be set within the program (0 by default), where a taxonomic classification is only taken if it is above this pre-defined threshold. For a short example of this, have a look at ([this figure in the Kraken2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46#Fig1)).

Here, we can see that the query sequence has been broken up into 27 k-mers. 11 of these k-mers are unclassified (shown in white), while 10 have been classified as the blue taxon, 4 as the orange taxon, 1 as the purple taxon and 1 as the red taxon. In this example, the purple taxon is the ancestor of all taxa and the blue is the ancestor of the orange. The sequence will be classified to the leaf node, or the lowest rank possible, and we can use the number of k-mers that were classified to calculate how confident we are in that taxonomic assignment. Our confidence that this sequence belongs to the orange taxon if 4/27, or 0.15. Our confidence that it belongs to the blue taxon is (10+4)/27 = 14/27 = 0.52, and our confidence that it belongs to the purple taxon is the sum of all classified k-mers, so 16/27 = 0.59. By default, the confidence threshold that Kraken uses is 0.00, which would mean that the sequence would be classified as the orange taxon. If we increased the confidence threshold to 0.20 then the sequence would be classified as the blue taxon, and if we increased it to 0.60 or higher then the sequence would be unclassified.

Again, you can see [our paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000949) if you want to see how this impacts the classifications given in more detail (and for more guidance on this), but as the confidence threshold is increased, fewer reads will be classified, but we can be more certain that the reads that are classified are classified correctly. To some extent, the "best" confidence threshold to use will depend on your research question and whether you are more concerned with ensuring that you identify all taxa that you possibly can (no false negatives/high recall), or that all of the taxa that you do identify are actually in the sample (no false positives/high precision). If you are running your own samples, then it may be worth taking a small (representative) subset and running them with a few different thresholds to see how the confidence threshold impacts this.

### Back to the Kraken2 output

Now we will take a look at some of the output files from Kraken2. First take a look at a file in the ```kraken2_kreport``` folder with ```head -30 kraken2_kreport/BB190_0.0_minikraken.kreport```. The first few lines should look something like this:
```{bash}
97.23	4823613	4823613	U	0	unclassified
  2.77	137575	93	R	1	root
  2.77	137343	1184	R1	131567	  cellular organisms
  2.68	133086	12486	D	2	    Bacteria
  1.29	63812	5448	P	1224	      Pseudomonadota
  0.79	39351	3632	C	28211	        Alphaproteobacteria
  0.59	29301	3457	O	356	          Hyphomicrobiales
  0.38	18619	2121	F	41294	            Nitrobacteraceae
  0.31	15439	5852	G	374	              Bradyrhizobium
  0.07	3672	3672	S	1437360	                Bradyrhizobium erythrophlei
  0.05	2456	543	G1	2631580	                unclassified Bradyrhizobium
  0.00	170	170	S	2782654	                  Bradyrhizobium sp. 186
  0.00	162	162	S	2782665	                  Bradyrhizobium sp. 200
  0.00	131	131	S	2782641	                  Bradyrhizobium sp. 170
```

The columns here are (from the [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual)):
| Column Number      | Data Type     |
| :------------- | :----------: |
| 1 | The percentage of sequences covered by the clade rooted at this taxon |
| 2 | The number of sequences covered by the clade rooted at this taxon |
| 3 | The number of sequences assigned directly to this taxon |
| 4 | A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies |
| 5 | NCBI taxonomic ID number |
| 6 | Indented scientific name |

You should notice first that we have a very large proportion of unclassified reads. In the example output above, 97.23% of sequences were unclassified while only 2.77% of sequences were classified. The output is hierarchical, so going down, we can see that 2.77% of sequences were from cellular organisms (so almost all of the classified sequences), and 2.68% were Bacterial sequences. If we look at the second and third columns for the Bacteria, we can see that 12,486 of the 4,823,613 sequences that were classified as Bacteria weren't classified to any rank lower than Bacteria (Domain rank).

The reason that we have so few reads classified here is because we've used the smallest Kraken2 database. We typically find that as the size of the database increases, so does the proportion of reads that can be classified because there's a higher chance that something similar to each of your reads exists in the database. Unfortunately the full RefSeq database is far too big to use in this workshop (it's ~800 GB), but that's what we've used to classify the samples for the data that you're using in the Phyloseq part of the workshop, below. Now take a look at the same sample that's been classified using the full database: ```head -30 ~/CourseData/MIC_data/metagenome_data/kraken2_kreport/BB190_0.0_RefSeqV214.kreport```

**Question 1: How many more reads are classified using the full database than with the mini database?**

Now take a look at a sample that was classified using a different confidence threshold with the full database, *e.g.*: ```head -30 ~/CourseData/MIC_data/metagenome_data/kraken2_kreport/BB190_0.1_RefSeqV214.kreport```

**Question 2: What do you notice about the number, or proportion, of reads that have been classified?**

Now look at the same sample but classified with the mini database and a higher confidence threshold, *e.g.*: ```head -30 kraken2_kreport/BB190_0.1_minikraken.kreport```

**Question 3: What do you notice now about the number of reads that have been classified? Does this seem like an appropriate confidence threshold to use?**

## Running Bracken

The raw output files created by Kraken2 should be corrected using Bracken to get more accurate estimations on the abundance of different taxa at various taxonomic levels. While Kraken2 classifies each of our reads to their best locations within a taxonomic tree, in some cases this will result in a read being placed higher up in the tree due to two species or sometimes even genera (or higher ranks) sharing the exact same sequence (k-mer) to which that read matched. Bracken can be used to solve this issue by taking these reads and mapping them back to the genomes they might belong to. By doing this it generates new abundance estimates that take into account the problem laid out above.

Once Kraken2 has finished running, we can run Bracken on these samples. Bracken doesn't require much memory and usually runs very quickly. 

First we will make a directory to output the results to:
```{bash}
conda activate bracken
mkdir bracken_out
```

Now we'll run Bracken:
```
parallel -j 2 'bracken -d ~/CourseData/MIC_data/metagenome_data/k2_pluspf_08gb_20230314/ -i {} -o bracken_out/{/.}.species.bracken -r 100 -l S -t 1' ::: kraken2_kreport/*minikraken.kreport
```

Here we are again using Parallel to run bracken on any files in the ```kraken2_kreport/``` folder that end with ```minikraken.kreport```. 

Here, we're giving Bracken several options:
- ```-d``` - The database folder that contains the Bracken files in addition to the Kraken2 database (the e.g. ```database100mers.kmer_distrib``` files)
- ```-i``` - The input kreport file from kraken2, denoted for parallel by {1} because it's the first option outside of the quoted bracken command after the ```:::```
- ```-r``` - The length of the reads used for mapping
- ```-l``` - The level at which to estimate taxon abundance (S indicates species)
- ```-t``` - The number of reads required for a taxon to be used in the abundance estimate

Note that we generally want to remove some low abundance taxa (*i.e.*, increasing the ```-t``` argument), but I personally tend to leave these in at this point and will remove these in my downstream analyses. The reason for this is that k-mer based strategies such as kraken2 are notorious for making a large number of low abundance false positives. Check-out this [paper](https://www.cell.com/cell/fulltext/S0092-8674(19)30775-5) or our [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000949) for more information.

The final step in generating our taxonomic profiles is collecting all the results from each sample into a single file. This can be done using a helper script that comes installed with bracken:
```
mkdir bracken_out_merged
combine_bracken_outputs.py --files bracken_out/*minikraken.species.bracken -o bracken_out_merged/merged_output_minikraken.species.bracken
```

You might have noticed above that we activated a Bracken environment, so now we'll deactivate that and reactivate our general taxonomic one:
```
conda deactivate
conda activate taxonomic
```
It's often good practice to install bioinformatic tools inside their own conda environments. In this case, installing Bracken installs an older version of Kraken2 that we don't want to use, so we installed it in its own environment so it didn't affect us running Kraken2.

Lets take a look at our resulting file with the ```less``` command. You will notice that the file has different columns corresponding to different information about a single taxon.

| Column Name      | Data Type     |
| :------------- | :----------: |
| name |  Taxonomic ID  |
| taxonomy_id | NCBI taxonomic ID |
| taxonomy lvl | Character representing the level of taxonomy |
| SAMPLE.species.bracken_num | Number of reads that aligned to this classifcation for this SAMPLE|
| SAMPLE.species.bracken_frac | Proportion of reads that aligned to this classification for this SAMPLE |
| SAMPLE2.species.bracken_num | Number of reads that aligned to this classifcation for this SAMPLE2 |
| SAMPLE2.species.bracken_frac | Proportion of reads that aligned to this classification for this SAMPLE2 |
| etc.    |     |

We can also copy across the bracken files for the RefSeqCompleteV214 database, a script that we'll use now and our metadata:
```
ln -s ~/CourseData/MIC_data/metagenome_data/bracken_out_merged/merged_output_RefSeq.species.bracken bracken_out_merged/
ln -s ~/CourseData/MIC_data/metagenome_data/combine_merged_files.py .
ln -s ~/CourseData/MIC_data/metagenome_data/Blueberry_metadata_metagenome.tsv .
```

Finally, we'll create a single output file that contains the information for both the classifications using MiniKraken that we've run, as well as the RefSeqCompleteV214 database that we'd run prior to the workshop:
```
python combine_merged_files.py
```
This script also adds the full taxonomy information (domain, phylum, class, order, family, genus) for each of the species that Bracken has identified. You can take a look at it using ```head bracken_out_merged/bracken_combined.csv``` if you like. Note that this uses a file that isn't usually an output from Kraken2/Bracken - I've put it together based on the NCBI taxonomy information (which can be found [here](https://www.ncbi.nlm.nih.gov/taxonomy) for a browsable version or [here](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) for a downloadable version. It's not necessary for analysis of our samples, but I find it useful to have information on which e.g. family a given species belongs to and be able to collapse the species at higher ranks. It's at ```~/CourseData/MIC_data/metagenome_data/tax_file_species.csv```.

## Visualising the results using Phyloseq in R

Now we're going to explore the results a bit in RStudio. We'll mainly be using [Phyloseq](https://joey711.github.io/phyloseq/), which is a popular package for microbiome data analysis as once you've got your data into the correct format for it, it has lots of really useful functions that you can use for simple visualisation and exploration of your data. The initial steps of getting your data into the correct format for it can be a bit fiddly, but it's really useful and easy to use once you've performed these. 

Launch RStudio server and type in your username and password.

First of all, we'll import all of the packages that we'll be using:
```
library(phyloseq)
library(vegan)
library(ggplot2)
library(randomcoloR)
library(gridExtra)
library(tidyr)
```

Next we'll read in the file that contains our metadata and we'll also set the row names to be the SampleID and then remove the first column (so we don't have the row names duplicated):
```
metadata = read.csv('workspace/metagenome_workshop/Blueberry_metadata_metagenome.tsv', sep='\t')
rownames(metadata) = metadata$SampleID
metadata = metadata[,-1]
```

Next we'll read in the bracken output. We'll then separate out the column in this file that contains taxonomy information from the columns containing information about the abundance of each species in our samples.
```
bracken = read.csv('workspace/metagenome_workshop/bracken_out_merged/bracken_combined.csv')
bracken = as.data.frame(bracken)
row.names(bracken) = bracken$name

taxonomy = bracken[, c(1,50)]
row.names(taxonomy) = row.names(bracken)

table_num = data.matrix(bracken[,2:49])
rownames(table_num) = bracken[,1]
bracken = as.matrix(table_num)

taxonomy <- separate(data = taxonomy, col = Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\;")
taxmat <- taxonomy[,-1]
taxonomy = tax_table(taxmat)
taxa_names(taxonomy) <- rownames(taxmat)
```

Now we'll convert the bracken output, metadata and taxonomy information to a phyloseq object:
```
physeq = phyloseq(otu_table(bracken, taxa_are_rows = TRUE), sample_data(metadata), taxonomy)
```
Now we have our data imported to the phyloseq object (named "physeq") and from here we'll just be using this object.

First we'll take a look at the number of sequences within each sample have been classified using the different database and confidence threshold options:
```
sums = as.data.frame(sample_sums(physeq))
sums[,'Database'] = sample_data(physeq)$Database
sums[,'Confidence.Threshold'] = sample_data(physeq)$Confidence.Threshold
colnames(sums) = c('Sum', 'Database', 'Confidence.Threshold')
sums$Confidence.Threshold <- as.factor(sums$Confidence.Threshold)
ggplot(sums, aes(x=Confidence.Threshold, y=Sum, fill=Database)) + geom_boxplot() + scale_y_log10() + labs(y= "Number of reads classified", x = "Confidence Threshold")
```

**Question 4: How many reads are classified with each of the database/confidence threshold combinations? Which would you choose based on this?**

Now we'll just reduce the phyloseq object to look at the full database and we'll also perform the normalisations that we will want to use:
```
db = "RefSeqV214"
conf = 0.1
physeq_red <- prune_samples(sample_data(physeq)$Database == db, physeq) #keep only samples where the database is "RefSeqV214"
physeq_red <- prune_samples(sample_data(physeq_red)$Confidence.Threshold == conf, physeq_red) #keep only samples where the confidence threshold is 0.1
physeq_rare <- rarefy_even_depth(physeq_red, sample.size = min(sample_sums(physeq_red)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE) #rarefy to the lowest sample depth
physeq_relabun  <- transform_sample_counts(physeq_red, function(x) (x / sum(x))*100) #convert to relative abundance
```
There are different opinions in the microbiome community about how samples should be normalised to account for biases in sequencing, but we're not going into that in detail today. If you're interested and want to know more, you can find several papers linked in the [microbiome for beginners](https://github.com/LangilleLab/microbiome_helper/wiki/Microbiome-for-beginners) page.

### Stacked bar plots

We will often use something like a stacked bar plot to get a quick overview of the data. Here we'll use the relative abundance phyloseq object to get a quick overview of how the community is composed. 

#### Domain 

First we'll look at the domain (Bacteria/Archaea/Eukaryotes/Viruses) level:
```
palette = distinctColorPalette(30)
rnk = "ta1"
ps.rank = tax_glom(physeq_relabun, taxrank=rnk, NArm=FALSE)

plot_bar(ps.rank, fill=rnk) + facet_wrap(c(~Description_1, ~Description_3), scales="free_x", nrow=1) + theme(legend.text=element_text(size=5), legend.key.size = unit(0.3, "cm")) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=palette)
```
The first line here defines the colour palette that we want to use (and that we want it to have 30 distinct colours), the second line groups our phyloseq object to the domain level (called "ta1" here) and then the third line is used for plotting the stacked barplot. You can go ahead and run just parts of this line to see how each part of the code changes the plotted graph, e.g. ```plot_bar(ps.domain, fill="ta1")```, ```plot_bar(ps.domain, fill="ta1") + facet_wrap(c(~Description_1, ~Description_3), scales="free_x", nrow=1)```, etc.

**Question 5: Here you should see that there are also Eukaryotes in your samples. What do you notice about the distribution of the Eukaryotes in different samples?**

#### Phylum

Now we'll look at the phylum level:
```
palette = distinctColorPalette(30)
rnk = "ta2"
ps.rank = tax_glom(physeq_relabun, taxrank=rnk, NArm=FALSE)
rank.sum = tapply(taxa_sums(ps.rank), tax_table(ps.rank)[, rnk], sum, na.rm=TRUE)
top30 = names(sort(rank.sum, TRUE))[1:30]
ps.rank = prune_taxa((tax_table(ps.rank)[, rnk] %in% top30), ps.rank)

plot_bar(ps.rank, fill=rnk) + facet_wrap(c(~Description_1, ~Description_3), scales="free_x", nrow=1) + theme(legend.text=element_text(size=5), legend.key.size = unit(0.3, "cm")) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=palette)
```
Here we're doing the same as above, but we're doing this at the rank of "ta3" rather than "ta1", and because we now have more than just three taxa at this rank, we're first filtering our table to include only the 30 most abundant taxa. So the line starting with rank.sum is first calculating the sum of each phylum across all samples, the top30 line is sorting these sums to find out which the top 30 most abundant phyla are and then the ps.phylum = prune_taxa line is filtering the ps.phylum object to include only those top 30 most abundant taxa before we plot them. 

**Question 6: Can you modify this phylum level code to work at a lower rank?**
*Hint: You will need to modify the ta2 part*

### Alpha diversity

Now we'll have a look at some measures of alpha diversity in our samples. We're going to be looking at Richness (observed taxa), Chao1 richness (estimated richness) and Simpson's diversity (a measure of diversity which takes into account the number of species present, as well as the relative abundance of each species). The diversity measures can't be calculated on numbers that aren't integers, so we can't calculate these for the relative abundance data. We also don't have a tree for the metagenomic data so we can't use Faith's phylogenetic diversity.

**Question 7: Why do we not usually have a phylogenetic tree for read-based metagenomic analyses?**

First we'll plot alpha diversity individually for each sample:
```
plot_richness(physeq_rare, measures=c("Observed", "Chao1", "Simpson", "Shannon"))
```
Here there are other alpha diversity metrics that we could look at, but we're just focusing on these 3 that give a good overview of what our samples are like. 

And now we'll take a look at the differences in alpha diversity between the groups.

Per group (bulk vs rhizosphere):
```
plot_richness(physeq_rare, x="Description_1", measures=c("Observed", "Chao1", "Simpson", "Shannon")) + geom_boxplot()
```

Per group (forest vs managed):
```
plot_richness(physeq_rare, x="Description_3", measures=c("Observed", "Chao1", "Simpson", "Shannon")) + geom_boxplot()
```

**Question 8: What patterns do you notice with the alpha diversity between groups?**

### Beta diversity ordination plots

As you heard earlier, there are a lot of different beta diversity metrics that we could use, but we'll just take a look at two of them; Bray-Curtis dissimilarity (which takes into account the relative abundance of taxa) and Jaccard distance (which accounts only for the presence/absence of taxa). Here we're also using the ```adonis2``` function to perform PERMANOVA tests. These tell us about any significant differences between groups.

Bray-Curtis dissimilarity:
```
ps = physeq_rare
ps.ord <- ordinate(ps, "PCoA", "bray")
plot_ordination(ps, ps.ord, type="samples", color="Description_1", shape="Description_3")
distance <- phyloseq::distance(ps, method="bray", weighted=F)
adonis2(distance ~ sample_data(ps)$Description_1*sample_data(ps)$Description_3)
```
Here the lines are: 1) creating a copy of the physeq_rare object, 2) creating an ordination of this object, telling the function that we want to use the PCoA ordination method with the bray (short for Bray-Curtis dissimilarity) method, 3) plotting this ordination, using Description_1  to colour them and Description_3 for the shape of the points, 4) calculating a distance matrix with the bray method, and 5) performing a statistical test using the adonis2 function (PERMANOVA test) using Description_1 and Description_3 in the metadata.

**Question 9: What can you tell about the samples from how they group on the PCoA plots? Does this agree with the results of the PERMANOVA tests?**

Jaccard distance:  

```
ps = physeq_rare
ps.ord <- ordinate(ps, "PCoA", "jaccard")
plot_ordination(ps, ps.ord, type="samples", color="Description_1", shape="Description_3")
distance <- phyloseq::distance(ps, method="jaccard", weighted=F)
adonis2(distance ~ sample_data(ps)$Description_1*sample_data(ps)$Description_3)
```
Here we're doing exactly the same as above, but with Jaccard distance rather than Bray-Curtis dissimilarity.

**Question 10: How do these results compare with those for Bray-Curtis dissimilarity? What does this tell you about the abundances of different taxa?**

## Functional annotation introduction

The main goal of this part of the tutorial is to introduce students to functional profiling of taxonomic reads using MMSeqs and to visualise the results of both the taxonomic and functional annotations. As in the taxonomic composition part of this workshop, we want to emphasize that there is not a one-size-fits-all pipeline for analyzing MGS data, and another commonly used functional profiling tool is [HUMAnN 3](https://huttenhower.sph.harvard.edu/humann/).

Throughout this tutorial there will be questions are aimed at helping students understand the various steps in reference based metagenomic profiling. [The answers are found on this page]().

## 1. Determining functional profiles using MMseqs2

Now that we have an idea of the taxonomic profiles of our samples we can also look at their functional potential. To do this we will being using [MMseqs2](https://github.com/soedinglab/mmseqs2) along with a few scripts developed by our lab to connect the functional annotations to the taxonomic annotations.

MMseqs2 works by taking sequenced reads, translating them into protein and then mapping them against a protein database. In this case we will be using [UniRef90](https://www.uniprot.org/help/uniref), a large protein database clustered at 90% identity. Our lab has developed a set of scripts to run this tool and use the information to assign both a taxonomic ID and a functional ID to each read within our samples.

As for the taxonomic composition, we would usually perform some initial steps of quality control, removal of unwanted reads and joining of forward and reverse reads, but for now, we've done this for you. 

We already created the directory that we'll be using in the taxonomic part of this workshop, so we'll change into that and as previously, we've installed all of the tools that we'll need for this part of the tutorial into an Anaconda environment, which we can activate with the following commands:

```
cd workspace/metagenome_workshop
conda activate functional
```

Now we'll create a folder for our MMSeqs output to go into:
```
mkdir mmseqs_U90_out
```

The first step will be running a command creates an MMSeqs database from the the input fastq files. The creation of this database is necessary for MMSeqs as it vastly increases the speed at which translated DNA sequences can be mapped against a protein database:
```
parallel -j 4 --progress 'mmseqs createdb {} mmseqs_U90_out/mmseqs-{/.}queryDB' ::: cat_reads/*
```
Notice that we're using Parallel for this again, and this time we're running 4 jobs at once and we've used the ```--progress``` flag to see which job we're on and how many have already been run. There aren't many options in the command, but you can see that after ```createdb``` this is where the reads file is taken in, and then we modify this name slightly for the output that will go into the ```mmseqs_U90_out``` folder.

The next commands we unfortunately won't be able to run as they require more memory than we have, and would take a very long time to run (I would usually assume >24 hours to run this for a set of samples), but we'll go through what they're doing still.

This command is the real meat of the job file and runs the freshly created sample database against the provided UniRef90 protien database:
```
parallel -j 1 'mmseqs search mmseqs_U90_out/mmseqs-{/.}queryDB ~/CourseData/MIC_data/metagenome_data/MMSeqs2_db/mmseqsUniref90DB mmseqs_U90_out/mmseqs-{/.}resultDB tmp --db-load-mode 3 --threads 4 --max-seqs 25 -s 1 -a -e 1e-5 > /dev/null 2>&1' ::: cat_reads/*
```
*Note: If you ran this command by accident, you can stop it again by pressing ```ctrl```+```c``` (at the same time)*

There are a number of parameters in this command:
* ```--db-load-mode 3``` - This parameter tells MMSeqs how to deal with loading the database into memory. For more information you can check out this [page](https://mmseqs.com/latest/userguide.pdf). However, setting this parameter to 3 helps when running MMSeqs on a cluster environment. 
* ```--threads``` - The number of processors we want MMSeqs to use during the search
* ```--max-seqs 25``` - This indicates that we want MMSeqs to output at maximum 25 hits for each sequence
* ```-s 1``` - This indicates the sensitivity that we want MMSeqs to run at. Increasing this number will lower the speed at which mmseqs runs but will increase its sensitivity. For well-explored environments such as the human gut, a setting of 1 should suffice.
* ```-a``` - This indicates that we want our results to output backtraces for each sequence match. These are needed to convert the resulting MMSeqs file into a usable file format. 
* ```-e 1e-5``` - This indicates that we only want to keep matches that are below an E-value of 1e-5 (E-values are a measure of how well two sequences match one another, and the closer they are to zero, the better the match is).
* ```> /dev/null 2>&1``` This final part of the command allows us to run the command without having too much text printed to our screen.

The final command allows us to convert the resulting file from the MMSeqs2 format into one that is more usable:
```
#parallel -j 1 'mmseqs convertalis mmseqs_U90_out/mmseqs-{/.}queryDB MMSeqs2_db/mmseqsUniref90DB mmseqs_U90_out/mmseqs-{/.}resultDB mmseqs_U90_out/mmseqs-{/.}-s1.m8 --db-load-mode 2 > /dev/null 2>&1' ::: cat_reads/*
```
This command is similar and takes as input the query database we made from our first command, the UniRef90 database we searched against and the resulting file from our search command. It will output the file ```mmseqs_U90_out/*.m8```. 

Now we'll get the results that would have been created from these steps.
We'll create a new folder for the ```.m8``` files and create a link to them:
```
mkdir mmseqs_m8_files
ln -s ~/CourseData/MIC_data/metagenome_data/mmseqs_U90_out/*.m8 mmseqs_m8_files/
```
The ```mmseqs search``` also made a lot of files that we're not interested in! We didn't copy them across here, but you can see them by running ```ls ~/CourseData/MIC_data/metagenome_data/mmseqs_U90_out/``` If you were running this for yourself, after you have the ```.m8``` files you could remove the others to save significantly on the amount of hard drive space that is being taken up by these files.

Lets take a quick look at one of the files we just moved into the directory ```mmseqs_m8_files``` using the less command.

```
less mmseqs_m8_files/mmseqs-BB209-s1.m8
```

We you will see is a file in BLAST tabular format.

| Column Number       | Data Type     |
| :------------- | :----------: |
| 0 |  query sequence ID  |
| 1 | Subject (database) sequence ID |
| 2 | 	Percent Identity |
| 3 | Alignment Length |
| 4 | Number of gaps |
| 5 | 	Number of mismatches |
| 6 | Start on the query sequence |
| 7 | End on the query sequence |
| 8 | Start on the database sequence |
| 9 | 	End on the database sequence |
| 10 | 	E value - the expectation that this alignment is random given the length of the sequence and length of the database |
| 11 | bit score - the score of the alignment itself |


**Question 1: How many protein sequences did the sequence ```SRR8742630.234641``` align with in the sample BB198? Which alignment/alignments have the lowest E-value/highest bitscore?**

You can press ```q``` to get out of the less view again.

The next step we need to take is to get the name of the protein sequence that had the best alignment for each sequence read in our samples. We can achieve this by running the command:
```
mkdir mmseqs_U90_out_tophit
cp -r ~/CourseData/MIC_data/metagenome_data/Functional_Helper_Scripts/ .
python Functional_Helper_Scripts/pick_uniref_top_hit.py --unirefm8Dir mmseqs_m8_files --output_path mmseqs_U90_out_tophit
```

Now we have the best protein sequence that matches best with each of the sequences in our samples. At this point, we could just look at the distribution of functions in our samples by adding the number of hits to each sequence within our samples, and we could analyse this in the same way that we looked at the taxonomic annotation previously. But we think it's nice to combine the taxonomic and the functional annotations together. 

So we'll be creating a final data table that contains the stratified abundance of each function in our samples. The next steps are using various databases to create files that contain information on the taxonomic and functional annotation that we've obtained for each read in our initial sequence files as well as collating that information so that we know about the abundance of different taxa/function combinations within samples.

First we'll copy across a file that contains mapping between the raw kraken2 output files (that contain information on a read-by-read basis) and the MMSeqs2 output files for each sample:
```
cp ~/CourseData/MIC_data/metagenome_data/multi-sample-outfiles-w-m8.txt .
```
Have a look at it with the ```head``` or ```less``` commands if you like. 

Now we'll create a link to some of the database files within out current directory:
```
ln -s ~/CourseData/MIC_data/metagenome_data/MMSeqs2_db/* .
```
These databases contain information about the length of each gene in our UniRef protein database, which is important to normalise the abundance of each functional classification. 

The next step would be to run a script that will take in the taxonomy and function for each read and convert this to a more usable format:
```
python Functional_Helper_Scripts/parse_TaxonomyFunction.py --multisample multi-sample-outfiles-w-m8.txt --outputf Workshop_strat_matrix_RPKM.txt --stratified Y --map2EC Y
```
You can run this command if you like, but it takes quite a while to run, so for the sake of time it is probably easier to just grab the file that we've already made:
```
ln -s ~/CourseData/MIC_data/metagenome_data/Workshop_strat_matrix_RPKM.txt .
```
This command would have generated this final stratified data matrix that shows the abundance of each EC (Enzyme Commission) number stratified by the different taxonomic classifications within the sample. This script also normalises the abundances of each of these ECs into reads per kilobase per million (RPKM). This abundance metric is the number of reads that mapped to that EC number per million reads mapped within the sample, divided by the gene length of that EC. It's important that the abundances are normalised by gene length or there would be an unfair bias toward longer genes due to the higher chance of them being sequenced. We can also run the same command as above without the ```--stratified Y``` option to generate functional profiles that are not broken down by the contributing taxa. 

Lets take a look at this file with the less command:
```
less Workshop_strat_matrix_RPKM.txt
```
*Remember that you can press ```q``` to exit the less view again*.

**Question 2: What is the RPKM contributed to the sample BB209 for the EC:5.4.2.2 contributed by Impatiens glandulifera?**

## 2. Adding descriptions and stratifying the file

Next, we want to add some descriptions to the EC numbers file and also convert this to pathway abundances. For this, we'll use some of the files included within [PICRUSt2](https://github.com/picrust/picrust2/wiki), but we've modified the scripts for this workshop.

First, we'll add the descriptions to our file:
```
ln -s ~/CourseData/MIC_data/metagenome_data/add_descriptions.py .
python add_descriptions.py
```
After running this, you should see a new file called ```Workshop_strat_matrix_RPKM_des.txt```. Take a look at it with the ```head``` command. You should see that each EC number now has a description associated.

**Question 3: What is the name of the enzyme with the EC number EC:5.4.2.2?**

The final step now is converting this into a stratified file and also adding the full taxonomy information to it. A stratified file is one where each line has a value for only the contribution of a gene by a single taxon in a single sample. Prior to this, we'll also take only the top most abundant pathways so that we'll be able to visualise it easily. If you were performing your own analysis, you might choose to visualise only *e.g.* those that were found to be significantly differentially abundant between your treatment and control samples. 


```
ln -s ~/CourseData/MIC_data/metagenome_data/add_taxonomy_and_stratify.py .
python add_taxonomy_and_stratify.py
```

## 3. Visualising with JARRVIS

If you've made it this far in the workshop - congrats!

If you didn't get through the rest of the workshop but still want to look at the end output, you can copy across the final file to your folder like so:
```
cp ~/CourseData/MIC_data/metagenome_data/EC-stratified-SankeyFormat.txt .
```

For this, the easiest thing to do is to download both of the files that we'll need. First we'll run:
```
cp ~/CourseData/MIC_data/16S_data/Blueberry_metadata_reduced.tsv .
sed -i 's/SampleID/sample_id/g' Blueberry_metadata_reduced.tsv
```
We've also just changed "SampleID" to "sample_id" as the first column name - often we find in bioinformatics tools they are very sensitive to small changes, and the tool we'll be using expects "sample_id" in this column rather than "SampleID". 

And then you can download these files from the ```metagenome_workshop/``` directory on http://##.uhn-hpc.ca/:
```
Blueberry_metadata_reduced.tsv
EC-stratified-SankeyFormat.txt
```

Next, we'll be using RStudio on your own computer. If you already have this installed, great. If not, download and install it from [here](https://posit.co/downloads/) (choose the free version download link!).

Go ahead and open it up and click into the console section. 

First, we're going to install a package called ```shiny```. You can do that with this command:
```
install.packages('shiny')
```
If you already had RStudio and it's already installed, running this won't do any harm.

Now we can run:
```
library(shiny)
runGist("943ff5fdbd94815cc27f302d9f56ff0b")
```
You should see a lot of different lines of code running, and eventually a window should pop up. This is a new tool that we're developing in the Langille Lab called [JarrVis](https://github.com/dhwanidesai/JarrVis) that we'll be using to visualise the links between the taxonomy and function in our different samples groups. 

You can now select the ```EC-stratified-SankeyFormat.txt``` file for the first box and the ```Blueberry_metadata_reduced.tsv``` file for the second box. Click ```Select Metadata Categories``` and it should populate the dropdown with the Metadata categories in our data. Choose either "Description_1" or "Description_3" here.

Leave ```Filtercollapse the dataframe``` as it is ("yes") and select a taxonomy level (I suggest starting with family or order). You can then click "Update the Gene Contribution Threshold from data". You can play around with the thresholds to use, but a lower end of around 100 gives a sensible number of things in the plot to look at. 

If you then click "Display Plot" then you should see a diagram come up that has sample labels on the left (red nodes), taxa names in the middle (blue nodes) and function names on the right (green nodes). Here, there are links between each of these - you can hover over things to highlight the full path and explore the different functions in your data and the taxa that contribute to each of these functions in your samples. 

