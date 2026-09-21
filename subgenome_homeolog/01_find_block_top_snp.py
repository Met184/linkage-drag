#!/usr/bin/env python3
"""01_find_block_top_snp.py

Most significant SNP per block, extended by 100 bp on each side.

Usage: python3 01_find_block_top_snp.py [--qtl FILE] [--snp FILE] [--out FILE]
"""

import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--qtl', default=os.path.join(HERE, 'QTL.txt'),
                    help='QTL block table (default: QTL.txt next to this script)')
parser.add_argument('--snp', default=os.path.join(HERE, '0SNP0.05N5.44.csv'),
                    help='significant SNP table (default: 0SNP0.05N5.44.csv next to this script)')
parser.add_argument('--out', default=os.path.join(HERE, 'QTL_with_topSNP.txt'),
                    help='output table (default: QTL_with_topSNP.txt next to this script)')
args = parser.parse_args()

QTL_FILE = args.qtl
SNP_FILE = args.snp
OUT_FILE = args.out

FLANK = 100  # extend 100 bp on each side

# ---------- 1. read SNPs, grouped by chromosome ----------
# Each SNP is kept as a dict for field access; p_wald is parsed as float for
# comparison.
snps_by_chr = {}
with open(SNP_FILE, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        chr_ = row['chr'].strip()
        try:
            ps = int(row['ps'])
            p = float(row['p_wald'])
        except (ValueError, KeyError):
            continue
        snps_by_chr.setdefault(chr_, []).append({
            'rs': row['rs'].strip(),
            'ps': ps,
            'p': p,
            'p_wald': row['p_wald'].strip(),
            'logP': row['logP'].strip(),
            'source': row['Source'].strip(),
        })

# ---------- 2. read QTL blocks ----------
blocks = []
with open(QTL_FILE, newline='', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    for row in reader:
        if len(row) < 6:
            continue
        qtl, block, chr_, start_s, end_s, trait = row[0], row[1], row[2], row[3], row[4], row[5]
        try:
            start = int(start_s)
            end = int(end_s)
        except ValueError:
            continue
        blocks.append({
            'qtl': qtl, 'block': block, 'chr': chr_,
            'start': start, 'end': end, 'trait': trait,
        })

# ---------- 3. most significant SNP per block ----------
# Output header: the original 6 columns plus 7 appended columns
out_header = ['QTL_IDl', 'Block_ID', 'chr', 'block_start', 'block_end', 'Trait',
              'topSNP_rs', 'topSNP_pos', 'topSNP_p', 'topSNP_logP', 'topSNP_source',
              'new_start', 'new_end']

n_no_snp = 0
out_rows = []
for b in blocks:
    candidates = []
    for s in snps_by_chr.get(b['chr'], []):
        if b['start'] <= s['ps'] <= b['end']:
            candidates.append(s)
    if not candidates:
        n_no_snp += 1
        out_rows.append([b['qtl'], b['block'], b['chr'], str(b['start']), str(b['end']), b['trait'],
                         '', '', '', '', '', '', ''])
        continue
    top = min(candidates, key=lambda s: s['p'])  # smallest p = most significant
    new_start = top['ps'] - FLANK
    new_end = top['ps'] + FLANK
    out_rows.append([b['qtl'], b['block'], b['chr'], str(b['start']), str(b['end']), b['trait'],
                     top['rs'], str(top['ps']), top['p_wald'], top['logP'], top['source'],
                     str(new_start), str(new_end)])

# ---------- 4. write result ----------
with open(OUT_FILE, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(out_header)
    writer.writerows(out_rows)

print(f"total QTL blocks: {len(blocks)}")
print(f"blocks with a significant SNP: {len(blocks) - n_no_snp}")
print(f"blocks without any SNP in the interval: {n_no_snp}")
print(f"written: {OUT_FILE}")
