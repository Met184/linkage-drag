# Linkage drag load

| Script | Purpose |
|--------|---------|
| `01_weighted_drag_load.R` | per-accession QTL-internal and QTL-between drag load |
| `02_weighted_R_load.R` | per-accession antagonistic (R) load |

Three drag classes:

- **R (antagonistic)**: the accession carries the repulsion haplotype at a QTL
  classified as high-antagonism or low-antagonism according to its mean R rate.
- **QTL-internal drag**: two blocks within one QTL carry opposite-effect
  haplotypes for two traits (S-I or I-S).
- **QTL-between drag**: the same configuration formed by two different QTLs.

Each drag pair is weighted by `|d1| x |d2|`, the product of the absolute effect
sizes of the two traits involved, where `d` is Cohen's d between the accessions
carrying the superior and the inferior haplotype of a QTL-trait combination.
Weights are normalized to mean 1, and only pairs with a pure-drag rate above 20%
are kept.

For every accession an unweighted count (each pair counts 1) and a weighted sum
are reported, both as absolute values and as rates over the number of valid pairs.

Input: `encode/*_SIR.csv`, `bio/<TRAIT>.csv`, `antagonism/antagonism_by_QTL.csv`,
`internal_drag/QTL_trait_pairs.csv`, `between_drag/QTL_between_drag_details.csv`
and `meta/yangbenALL.csv`.

Requires R with `data.table`, `ggplot2`, `rstatix`, `multcompView`.

```bash
Rscript 01_weighted_drag_load.R --work DIR --out DIR
Rscript 02_weighted_R_load.R    --work DIR --out DIR
```
