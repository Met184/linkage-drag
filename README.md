# Linkage drag and genomic compensation in global cotton improvement

Analysis code for the manuscript "Breeding-Driven Dynamics of Linkage Drag and
Genomic Compensation in Global Cotton Improvement".

```
haplotype_coding/      block-wise haplotype coding
epistasis/             epistasis between the A07 drag QTLs and remote FL/FS QTLs
subgenome_homeolog/    subgenome homeolog analysis of genomic compensation
drag_load/             load of the three classes of linkage drag
recombination/         relative recombination rate across three populations
```

Each folder has a short `README.md`. Every script takes its paths from the command
line, so no absolute path needs to be edited.

Population labels follow the manuscript (Group-A ... Group-G), defined in the Data
collection section of the manuscript. They are written `Group-B` in data values
and file names and `Group_B` in R identifiers.

Raw genotype data, phenotype tables and result tables are not included.
