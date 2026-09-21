# Subgenome homeolog analysis (genomic compensation)

Scripts implementing the subgenome homeolog analysis that separates the two
components of genomic compensation:

- **subgenome buffering**: the homeologous copy of a drag block carries a
  superior haplotype in the other subgenome and therefore provides structural
  redundancy for the locally inferior haplotype;
- **polygenic compensation**: the homeologous interval is itself inferior, so
  the phenotypic penalty of the drag block is offset by superior alleles at
  other, unlinked loci.

## Pipeline

| Step | Script | Purpose | Runs on |
|------|--------|---------|---------|
| 01 | `01_find_block_top_snp.py` | most significant SNP per block, +/- 100 bp -> 201 bp query intervals (`QTL_with_topSNP.txt`) | local |
| 02 | `02_blast_cross_subgenome.sh` | cross-subgenome BLAST (A segments against the whole D subgenome and vice versa) plus the reciprocal BLAST back to the original subgenome | server with BLAST+ and samtools |
| 03 | `03_parse_homeolog_rbh.py` | parse the forward and reverse BLAST output, call reciprocal best hits -> `homeolog_rbh_results.txt` | local |
| 04 | `04_count_block_snps.sh` | count variants and SNPs per block VCF -> `snp_count.tsv` | server |
| 05 | `05_homeolog_haplotype.R` | haplotype-phenotype association of the homeologous intervals (geneHapR, t-tests, coding) | local |
| 06 | `06_compensation_2x2.R` | 2x2 contingency table between the original block and its homeologous interval, with Fisher's exact test and the odds ratio | local |
| 07 | `07_compensation_by_group.sh` | the same contingency table for all accessions, for cultivated and for semi-wild cotton separately | local |
| 08 | `08_sequence_similarity.py` | sequence similarity of the homeologous intervals, block level and per sample | server (samtools, bcftools) |
| 09 | `09_effect_size.R` | effect size of the original block and of its homeologous interval, and the ratio between the two | local |

## Definition of the compensation call

For every significant block the 2x2 table crosses the coding of the original
block with the coding of its homeologous interval:

    SS = original superior & homeolog superior
    SI = original superior & homeolog inferior
    IS = original inferior & homeolog superior
    II = original inferior & homeolog inferior

    pct_homeo2_given_orig2 = SS / (SS + SI)

    > 0.5  -> concordant_homeolog_superior   (subgenome buffering)
    < 0.5  -> compensated_homeolog_inferior  (polygenic compensation)
    = 0.5  -> no_preference_0.5

Coding rule (identical to the main haplotype coding script): 2 = superior,
0 = inferior, and each haplotype is coded relative to the mean of the haplotype
phenotypic means within the block. `VW` and `FD` are reverse-coded because
higher values are agronomically undesirable.

Odds ratios are reported for reference: OR > 1 indicates a positive
(concordant) association and OR < 1 a negative (compensated) association.

## Inputs

| Input | Used by |
|-------|---------|
| `QTL.txt`, `0SNP0.05N5.44.csv` | 01 |
| reference genome FASTA + `.fai` | 02, 08 |
| `QTL_with_topSNP.txt` | 02, 03, 08 |
| `homeolog_rbh_results.txt` | 05, 08 |
| `block/Block_XXX_L.vcf` | 04, 05, 06, 07 |
| `consensus_coding/QTL_XXX/{trait}_consensus_coding.csv` | 06, 07 |
| `phenotypes/{trait}.csv` (columns `SRR`, `bio`, `<TRAIT>`) | 05, 06, 07 |
| `haplotype_means_coding.csv`, `haplotype_ttest_results.csv`, `haplotype_sample_code.csv`, `Block_XXX.vcf` | 08 |
| master VCF with all samples | 08 |
| per-accession detail tables of the selected blocks | 09 |

## Requirements

- Python 3 (standard library only; `samtools` / `bcftools` on `PATH` for step 08)
- R (>= 4.5.0) with `geneHapR`, `dplyr`, `tidyr`, `stringr` for steps 05, 06 and 09
- BLAST+ (2.17.0) and samtools (1.17) for steps 02, 04 and 08

## Configuration

Every script takes its paths from the command line, or from environment variables
for the shell scripts, so no absolute path needs to be edited.

```bash
python3 01_find_block_top_snp.py --qtl QTL.txt --snp 0SNP0.05N5.44.csv

GENOME=/path/genome.fasta WORK=/path/project bash 02_blast_cross_subgenome.sh

python3 03_parse_homeolog_rbh.py --qtl QTL_with_topSNP.txt

bash 04_count_block_snps.sh /path/project/block /path/project/snp_count.tsv

Rscript 05_homeolog_haplotype.R --home /path/project --ref /path/coding_project \
                                --vcf /path/project/block --out /path/project/haplotype_analysis

Rscript 06_compensation_2x2.R --home /path/project \
                              --consensus /path/coding_project/consensus_coding \
                              --ref /path/coding_project --out /path/project/compensation

BASE=/path/project bash 07_compensation_by_group.sh

python3 08_sequence_similarity.py --work /path/project \
        --genome /path/genome.fasta --master-vcf /path/all_samples.vcf.gz

Rscript 09_effect_size.R --home /path/project --detail /path/project/detail \
                         --out /path/project/effect_size
```

## Notes

- Population labels follow the manuscript (Group-A ... Group-G) and are used in
  the `bio` column of the phenotype tables, in the sample metadata and in the
  population column of the `COMBOS` definition.
- The block x trait combinations analysed in steps 06 and 07 are fixed in the
  `COMBOS` definition of those scripts, and the cases analysed in step 09 are
  fixed in the `cases` definition of that script.
- Raw genotype data, phenotype tables and all result tables are not part of this
  code repository. See the Availability of Data and Materials section of the
  manuscript.
- Figure scripts for this module are not included here.
