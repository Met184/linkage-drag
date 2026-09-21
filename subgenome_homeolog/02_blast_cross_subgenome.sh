#!/bin/bash
# 02_blast_cross_subgenome.sh
# Cross-subgenome BLAST of the block query intervals plus the reciprocal BLAST.
# Usage: GENOME=genome.fasta WORK=dir bash 02_blast_cross_subgenome.sh

set -e

# Cluster environment; skipped when `module` is not available
module load samtools/1.17 blast/2.17.0 2>/dev/null || true

GENOME="${GENOME:-/path/to/Ghirsutumv1.1_genome.fasta}"
WORK="${WORK:-$PWD}"
cd "$WORK"

# ---------- 1. split the A and D subgenome sequences ----------
grep -P '^Ghir_A' "$GENOME".fai | cut -f1 > A_chrs.txt
grep -P '^Ghir_D' "$GENOME".fai | cut -f1 > D_chrs.txt
samtools faidx "$GENOME" $(cat A_chrs.txt) > A_subgenome.fa
samtools faidx "$GENOME" $(cat D_chrs.txt) > D_subgenome.fa
echo "subgenome fasta: A=$(grep -c '>' A_subgenome.fa)  D=$(grep -c '>' D_subgenome.fa)"

# ---------- 2. build BLAST databases ----------
makeblastdb -in A_subgenome.fa -dbtype nucl -out A_db -title A_db >/dev/null
makeblastdb -in D_subgenome.fa -dbtype nucl -out D_db -title D_db >/dev/null

# ---------- 3. extract queries, split by subgenome ----------
> query_A.fa
> query_D.fa
tail -n +2 QTL_with_topSNP.txt | tr -d '\r' | while IFS=$'\t' read -r qtl blk chr bs be trait rs pos p logP src ns ne; do
  out="query_A.fa"; case "$chr" in Ghir_D*) out="query_D.fa";; esac
  samtools faidx "$GENOME" "${chr}:${ns}-${ne}" | sed "1s/.*/>${blk}|${chr}/" >> "$out"
done
echo "queries: A=$(grep -c '>' query_A.fa)  D=$(grep -c '>' query_D.fa)"

# Helper: keep only the highest-bitscore hit per query
# (group by qseqid, sort by column 12 = bitscore descending, take the first)
top1() {
  sort -k1,1 -k12,12nr "$1" | awk '!seen[$1]++'
}

# ---------- 4. forward cross-subgenome BLAST (A->D, D->A), single-threaded ----------
blastn -query query_A.fa -db D_db -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 > fwd_A_to_D.raw.tsv
blastn -query query_D.fa -db A_db -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 > fwd_D_to_A.raw.tsv
top1 fwd_A_to_D.raw.tsv > fwd_A_to_D.tsv
top1 fwd_D_to_A.raw.tsv > fwd_D_to_A.tsv
echo "forward hits: A->D $(wc -l < fwd_A_to_D.tsv)  D->A $(wc -l < fwd_D_to_A.tsv)"

# ---------- 5. reverse queries (forward best hits only), blasted back to the original subgenome ----------
> rev_query_D.fa
awk -F'\t' '{s=$9; e=$10; if(s>e){t=s;s=e;e=t} print $2":"s"-"e"\t"$1}' fwd_A_to_D.tsv \
  | while IFS=$'\t' read -r region name; do
      samtools faidx D_subgenome.fa "$region" | sed "1s/.*/>${name}/"
    done > rev_query_D.fa
blastn -query rev_query_D.fa -db A_db -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 > rev_D_to_A.raw.tsv
top1 rev_D_to_A.raw.tsv > rev_D_to_A.tsv

> rev_query_A.fa
awk -F'\t' '{s=$9; e=$10; if(s>e){t=s;s=e;e=t} print $2":"s"-"e"\t"$1}' fwd_D_to_A.tsv \
  | while IFS=$'\t' read -r region name; do
      samtools faidx A_subgenome.fa "$region" | sed "1s/.*/>${name}/"
    done > rev_query_A.fa
blastn -query rev_query_A.fa -db D_db -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 > rev_A_to_D.raw.tsv
top1 rev_A_to_D.raw.tsv > rev_A_to_D.tsv

echo "reverse hits: D->A $(wc -l < rev_D_to_A.tsv)  A->D $(wc -l < rev_A_to_D.tsv)"
echo "done; next: python3 03_parse_homeolog_rbh.py"
