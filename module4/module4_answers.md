# Metagenome taxonomic annotation with Kraken2 and Phyloseq

**Question 1: How many more reads are classified using the full database than with the mini database?**

Using the Mini database, 137,575 reads were classified (2.77%), while using the full database, 1,974,579 reads were classified (39.80%), so ~1.8M more reads are classified using the full database.

**Question 2: What do you notice about the number, or proportion, of reads that have been classified?**

With a confidence threshold of 0.1 (using the full database), only 717,094 reads have been classified (14.45%), so as the confidence threshold increases, fewer reads are classified.

**Question 3: What do you notice now about the number of reads that have been classified? Does this seem like an appropriate confidence threshold to use?**

With a confidence threshold of 0.1 (using the mini database), only 2,770 reads have been classified (0.06%), so with this database too, increasing the confidence threshold leads to fewer reads being classified. This seems like too low a proportion of reads that are classified to get an accurate overview of the community, though, so if you see this with your own data, I'd suggest using a lower confidence threshold or a bigger database. The right database/confidence threshold combination to use will depend on your data as well as the research question at hand, and I'd often suggest running a few different confidence thresholds with a subset of your data to see what seems sensible.

**Question 4: How many reads are classified with each of the database/confidence threshold combinations? Which would you choose based on this?**

We can see that at a confidence threshold of 0, an average of approximately 10<sup>5</sup> or above 10<sup>6</sup> reads have been classified with the MiniKraken or RefSeqV214 databases, respectively, while at a confidence threshold of 0.1, an average of below 10<sup>4</sup> or just below 10<sup>6</sup> reads have been classified with the MiniKraken or RefSeqV214 databases, respectively. It seems that with using the MiniKraken database, the decrease in the number of reads classified is several orders of magnitude when moving from a confidence threshold of 0 to 0.1, which doesn't seem acceptable even if we have an increase in precision - if we have to use the MiniKraken database because that's all we have the resources for, we should use a confidence threshold of 0. For RefSeqV214, there is a decrease in the number of reads classified, but it's not a full order of magnitude and seems more reasonable to have an increase in precision, so I would choose the confidence threshold of 0.1. 

**Question 5: Here you should see that there are also Eukaryotes in your samples. What do you notice about the distribution of the Eukaryotes in different samples?**

In the bulk samples Eukaryotes make up about 10% of the community, while in the rhizosphere samples they make up 35-70% of the samples. They don't seem to vary much between forest and managed samples, although they are particularly high in a single managed rhizosphere sample.

**Question 6: Can you modify this phylum level code to work at a lower rank?**

An example of this at the family level:
```
palette = distinctColorPalette(30)
rnk = "ta5"
ps.rank = tax_glom(physeq_relabun, taxrank=rnk, NArm=FALSE)
rank.sum = tapply(taxa_sums(ps.rank), tax_table(ps.rank)[, rnk], sum, na.rm=TRUE)
top30 = names(sort(rank.sum, TRUE))[1:30]
ps.rank = prune_taxa((tax_table(ps.rank)[, rnk] %in% top30), ps.rank)

plot_bar(ps.rank, fill=rnk) + facet_wrap(c(~Description_1, ~Description_3), scales="free_x", nrow=1) + theme(legend.text=element_text(size=5), legend.key.size = unit(0.3, "cm")) + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=palette)
```
Note that here we've changed the ta2 (phylum) to ta5 (family; domain=ta1, phylum=ta2, class=ta3, order=ta4, family=ta5, genus=ta6, species=ta7).

**Question 7: Why do we not usually have a phylogenetic tree for read-based metagenomic analyses?**

In amplicon analyses, we're looking at reads from only a single gene, so we can make trees directly with these genes. In metagenomic analyses, the reads come from different parts of the genome and it would therefore be much more computationally intensive to try to build a tree with these. Often when we assemble the reads into contigs and then use these to build metagenome assembled genomes (MAGs) then we will perform tree-based analyses using this. However, there are some methods that help to generate trees for metagenomic data, e.g. [PhyloT](https://phylot.biobyte.de/index.cgi) - which will generate a taxonomic tree based on a list of NCBI taxonomy ID's - or mapping of the prokaryotic reads to genomes that have already been placed in a phylogenetic tree, such as that provided with each [Genome Taxonomy Database release](https://gtdb.ecogenomic.org/).

**Question 8: What patterns do you notice with the alpha diversity between groups?**

Richness (observed taxa and Chao1) is higher in rhizosphere than bulk samples, but Simpson's diversity is lower. Richness is slightly lower in managed than forest samples but Simpson's diversity is similar. This suggests that there are more taxa present in rhizosphere than bulk samples (higher richness), but that just a few of these taxa are dominant in the rhizosphere samples (lower diversity).

**Question 9: What can you tell about the samples from how they group on the PCoA plots? Does this agree with the results of the PERMANOVA tests?**

The samples group more strongly by whether they are from bulk or rhizosphere (Description_1) than whether they are from forest/managed (Description_3). We also see that the bulk/rhizosphere samples are separates on Axis 1, which accounts for a large amount of the overall variation between samples (81%).

You should see an output similar to this for the PERMANOVA test:
```
adonis2(formula = distance ~ sample_data(ps)$Description_1 * sample_data(ps)$Description_3)
                                                            Df SumOfSqs      R2       F Pr(>F)   
sample_data(ps)$Description_1                                1  0.66799 0.72628 24.8658  0.002 **
sample_data(ps)$Description_3                                1  0.02319 0.02522  0.8633  0.404   
sample_data(ps)$Description_1:sample_data(ps)$Description_3  1  0.01365 0.01484  0.5080  0.567   
Residual                                                     8  0.21491 0.23366                  
Total                                                       11  0.91973 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Here we can see the R<sup>2</sup> (effect size) and p (Pr(>F)) values for Description_1 (bulk vs rhizosphere), Description_3 (forest vs managed) and the interaction between the two (sample_data(ps)$Description_1:sample_data(ps)$Description_3). The residual shows us how much of the variation between samples is not explained by the variables that we're testing (23% in this case). So the stronger grouping on the PCoA plot by bulk vs rhizosphere is supported by the PERMANOVA test, where Description_1 accounts for 73% of the between sample variation (R<sup>2</sup>=0.72628) and is the only factor that is significant (*p*=0.002). Description_3 explains only a small amount of the variation between samples (2.5%) and this effect is not significant (p=0.404) and the interaction between Description_1 and Description_3 is also not significant (*p*=0.567) and accounts for only a small amount of the variation (1.5%).

**Question 10: How do these results compare with those for Bray-Curtis dissimilarity? What does this tell you about the abundances of different taxa?**

A lower proportion of the variation is explained by Jaccard distance than Bray-Curtis dissimilarity (Residual R<sup>2</sup>=0.32533), but the grouping of samples and results of the PERMANOVA test are similar. This suggests that the differences between the sample groups are driven more by presence/absence of taxa than they are by the abundance of those taxa.

# Metagenome functional annotation with MMSeqs2 and JarrVis

**Question 1: How many protein sequences did the sequence ```SRR8742630.234641``` align with in the sample BB198? Which alignment/alignments have the lowest E-value/highest bitscore?**

To have a look at this, we can use the ```grep``` command, like so: ```grep -c "SRR8742630.234641" mmseqs_m8_files/mmseqs-BB198-s1.m8``` - in this case, we've asked it to count how many times "SRR8742630.234641" appears in the file ```mmseqs_m8_files/mmseqs-BB198-s1.m8```. We should see that it appears 34 times. If we want to look at the alignments, we can use it again but without the ```-c``` (count) flag, so that it prints all of the alignments for us: ```grep "SRR8742630.234641" mmseqs_m8_files/mmseqs-BB198-s1.m8```. We should see that the lowest (best) E-value is the first match that comes up, and is 6.330E-18 (next to last column), and that this match also has the highest bitscore (last column). *Note: if you don't see 12 columns, you may need to reduce the text size in your terminal window/make the window wider*. 

**Question 2: What is the RPKM contributed to the sample BB209 for the EC:5.4.2.2 contributed by Impatiens glandulifera?**

You should be looking at the first column of the row named "EC:5.4.2.2|Impatiens glandulifera (taxid 253017)", so this is 8.813512568291685. *Note: remember that you can download these files easily and open them up in Excel (or another similar program)!*

**Question 3: What is the name of the enzyme with the EC number EC:5.4.2.2?**

Phosphoglucomutase (alpha-D-glucose-1,6-bisphosphate-dependent)
