# Haplotype coding

| Step | Script | Purpose |
|------|--------|---------|
| 01 | `01_encode_all_qtls.R` | code each LD-block haplotype as superior (2) or inferior (0) for every QTL, trait and population |
| 02 | `02_block_significance_by_population.R` | summarise significance and effect direction of every QTL x block x trait across populations and apply the block retention rule |

Step 01 produces the coding matrices from which Type I, Type II and Type III
linkage drag are derived. Coding rule: haplotypes below 1% frequency within a
population are discarded; each haplotype is compared with the mean of the
haplotype phenotypic means within the block (`>= mean` -> code 2, otherwise 0),
with the coding reversed for `VW` and `FD`, traits for which higher values are
undesirable.

Step 02 reads the step 01 output and reports, per QTL x block x trait, the
significance level and the effect direction in each population, the haplotypes
that are consistently superior or inferior, and the number of significant
populations. Its output is the table used for the retention rule applied in the
manuscript: blocks mapped in several populations are kept only when they are
significant in at least two populations with a consistent effect direction.

## Inputs

- step 01: `qtl_block_trait_map.csv` (`QTL_id`, `Block_ID`, `Trait_Code`),
  `vcf/Block_<id>.vcf`, `phenotypes/<TRAIT>.csv` (`SRR`, `bio`, `<TRAIT>`)
- step 02: the `recode_results/QTL_<id>/<TRAIT>/` directories written by step 01

The `bio` column carries the manuscript population labels (Group-A ... Group-G).

## Outputs

- step 01: `recode_results/QTL_<id>/<TRAIT>/` with the coding matrix, the
  haplotype phenotype means, the pairwise t-tests and a per-block significance
  summary
- step 02: `block_significance_summary.csv`

## Requirements

R with `geneHapR` (step 01 only), `dplyr`, `tidyr`, `stringr`.

```bash
Rscript 01_encode_all_qtls.R --base DIR --map FILE --vcf DIR --pheno DIR --out DIR
Rscript 02_block_significance_by_population.R --coding DIR --out DIR
```
