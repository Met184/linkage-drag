# Cross-QTL epistasis analysis

`cross_qtl_epistasis_analysis.R`

Tests epistatic interactions between the A07 drag QTLs (QTL051, QTL053, QTL054)
and genome-wide QTLs associated with fibre length (FL) and fibre strength (FS):

1. the A07 trans-eQTL catalogue is used to find the FL/FS QTLs on other
   chromosomes targeted by trans-eQTLs of genes in the QTL051-QTL054 intervals;
2. haplotype coding combinations are built for the significant QTL pairs;
3. the combinations are validated phenotypically and tested for epistasis.

Only accessions carrying either SH or IH at both loci are retained (n >= 30, and
n >= 20 with non-missing phenotypes). A two-way ANOVA `y ~ G_A * G_B` is fitted,
where `G_A` and `G_B` are the haplotype categories of the A07 QTL and of the
interacting QTL; the P-value of the interaction term is the test of epistasis.

Input: `qtl_block_significance.csv`, `eqtl.csv` and the per-QTL consensus coding
matrices `consensus_coding/QTL_<id>/<TRAIT>_consensus_coding.csv`.

Output: `trans_eQTL_FL_FS_summary.csv`, `cross_QTL_epistasis_summary.csv`,
`cross_QTL_combo_pheno.csv`, `cross_QTL_all_merged.csv`.

Requires R with `data.table` and `dplyr`.

```bash
Rscript cross_qtl_epistasis_analysis.R --work DIR --consensus DIR --out DIR
```
