#!/bin/bash
# 04_count_block_snps.sh
# Count the variants and SNPs of every Block_*_L.vcf in a directory.
# Usage: bash 04_count_block_snps.sh [VCF_DIR] [OUT_FILE]

set -e

VCF_DIR="${1:-.}"
OUT="${2:-snp_count.tsv}"
printf "Block\tTotal_variants\tSNPs\tOther\n" > "$OUT"

for f in "$VCF_DIR"/Block_*_L.vcf; do
  [ -e "$f" ] || continue
  blk="$(basename "${f%_L.vcf}")"
  # Single awk pass: total = number of non-comment lines;
  # snp = REF and ALT both single bases (simple biallelic SNP)
  read -r total snp < <(awk -F'\t' '
    !/^#/ {
      total++;
      if (length($4)==1 && $5 ~ /^[ACGTNacgtn]$/) snp++;
    }
    END { print total+0, snp+0 }
  ' "$f")
  printf "%s\t%d\t%d\t%d\n" "$blk" "$total" "$snp" "$((total-snp))" >> "$OUT"
done

nblock=$(( $(wc -l < "$OUT") - 1 ))
snpsum=$(awk -F'\t' 'NR>1{s+=$3} END{print s+0}' "$OUT")
echo "done: ${nblock} blocks, ${snpsum} SNPs in total"
echo "written: $OUT"
