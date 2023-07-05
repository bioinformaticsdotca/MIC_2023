## Microbial diversity metrics and data visualization with QIIME2

**Question 1: What is a good rarefaction depth for diversity analysis?**

A cut-off of 4,000 reads will be sufficient: the curve plateaus around this depth and we won't exclude any samples.

**Question 2: are there any significant differences in alpha diversity between any of our metadata categories?**

There are significant differences in richness and phylogenetic diversity between forest and managed environment samples. There are no significant differences between bulk and rhizosphere soil. 

**Question 3: which metadata category appears to provide more separation in the beta diversity PCoA plots?**

This is hard to say just by looking at the Emperor plots provided by QIIME2, but the forest/managed category appears to exhibit more distinct separation. You could verify this by exploring the `beta-group-significance` command within the `diversity` plugin.

**Question 4: can you identify any patterns between the metadata groups?**

Because stacked barcharts are limited in their analytical capabilities, it is hard to discern anything except very obvious patterns.

**Question 5: Does ANCOM identify any differentially abundant taxa between any of the metadata groups? If so, which one(s)?**

ANCOM does not identify any taxa as differentially significant between bulk and rhizosphere soils, but three ASVs are identified as differentially abundant between forest and managed environments. You can look up the identity of each ASV in the taxonomy.tsv file. 
