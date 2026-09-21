# Recombination rate

| Script | Purpose |
|--------|---------|
| `01_relative_recombination_rate.R` | per-QTL absolute and relative recombination rate in three populations, and the comparison of drag versus non-drag QTLs |
| `02_cross_population_correlation.R` | pairwise Spearman correlation of the relative rates on the common 100-kb windows |

Population-scaled recombination rates (rho = 4*Ne*r) depend on the SNP calling
pipeline and on population composition, so the absolute values differ between
populations. Cross-population comparison is therefore made on the **relative
recombination rate**:

    Rel = (rate of the interval or window) / (genome-wide mean rate)

`Rel = 1` equals the genome-wide mean; `Rel < 1` marks a cold region.

A QTL interval rate is the overlap-length weighted average of the windows it
spans. Genome-wide windows are obtained by aggregating each map into fixed bins
with the same weighting, keeping only bins present in all three populations.

Populations are labelled as in the manuscript: **Group-B, Group-C, Group-D**
(written `Group_B` etc. in R identifiers).

## Input maps

| File | What it is |
|------|------------|
| `maps/Group_B_fastEPRR_rho_50kb.txt` | raw FastEPRR output for Group-B, with header `chromosome Start End Rho CIL CIR`; `Rho` is the rate per 50-kb window and `CIL`/`CIR` its confidence bounds |
| `maps/Group_C_fastEPRR_rho_50kb.txt` | the same for Group-C |
| `maps/Group_D_published_cM_per_Mb.txt` | published Group-D map, without header, 4 columns `Chr Start End Recomb_Rate`, already in cM/Mb |
| `maps/Group_{B,C,D}_recombination_cM_per_Mb.csv` | the three maps rescaled to cM/Mb; the Group-B and Group-C rates are `Rho / 50` |

Script 01 also needs `QTL_intervals.csv` and the drag classification tables under
`drag/`.

Requires R with `data.table`.

```bash
Rscript 01_relative_recombination_rate.R --work DIR --out DIR
Rscript 02_cross_population_correlation.R --work DIR
```
