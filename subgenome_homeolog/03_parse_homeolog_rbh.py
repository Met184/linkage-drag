#!/usr/bin/env python3
"""03_parse_homeolog_rbh.py

Parse the forward and reverse BLAST output and call reciprocal best hits.

Usage: python3 03_parse_homeolog_rbh.py [--qtl FILE] [--out FILE] [--margin INT]
"""

import argparse
import os

HERE = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--qtl', default=os.path.join(HERE, 'QTL_with_topSNP.txt'),
                    help='QTL table with top SNPs (default: next to this script)')
parser.add_argument('--out', default=os.path.join(HERE, 'homeolog_rbh_results.txt'),
                    help='output table (default: next to this script)')
parser.add_argument('--margin', type=int, default=1000,
                    help='overlap tolerance in bp when calling an RBH (default: 1000)')
args = parser.parse_args()

QTL_FILE = args.qtl
OUT_FILE = args.out
MARGIN = args.margin
FWD_FILES = ['fwd_A_to_D.tsv', 'fwd_D_to_A.tsv']
REV_FILES = ['rev_D_to_A.tsv', 'rev_A_to_D.tsv']
DIR = os.path.dirname(os.path.abspath(QTL_FILE))

def strip(s):
    return s.strip().strip('\r')

# ---------- read the block table ----------
blocks = {}   # Block_ID -> dict
with open(QTL_FILE, encoding='utf-8') as f:
    f.readline()
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 13:
            continue
        blocks[strip(p[1])] = {'chr': strip(p[2]), 'trait': strip(p[5]),
                               'start': int(p[11]), 'end': int(p[12])}

def load_blast(files):
    """Return {qseqid (Block_ID): hit dict}; keep the highest bitscore."""
    d = {}
    for fn in files:
        path = fn if os.path.isabs(fn) else os.path.join(DIR, fn)
        if not os.path.exists(path):
            continue
        with open(path, encoding='utf-8') as f:
            for line in f:
                c = line.rstrip('\n').split('\t')
                if len(c) < 12:
                    continue
                qid = strip(c[0])
                blk = qid.split('|')[0]
                hit = {'sseqid': strip(c[1]), 'pident': float(c[2]),
                       'length': int(c[3]), 'sstart': int(c[8]), 'send': int(c[9]),
                       'evalue': float(c[10]), 'bitscore': float(c[11])}
                if blk not in d or hit['bitscore'] > d[blk]['bitscore']:
                    d[blk] = hit
    return d

fwd = load_blast(FWD_FILES)
rev = load_blast(REV_FILES)

out_header = ['Block_ID', 'Trait', 'query_chr', 'query_start', 'query_end',
              'homeolog_chr', 'homeolog_start', 'homeolog_end', 'strand',
              'pident', 'evalue', 'bitscore', 'RBH', 'rev_chr', 'rev_start', 'rev_end']
rows = []
n_fwd = n_rbh = 0
for blk, info in blocks.items():
    qchr = info['chr']
    qs, qe = info['start'], info['end']
    f = fwd.get(blk)
    if f is None:
        rows.append([blk, info['trait'], qchr, qs, qe, '', '', '', '', '', '', '', '', '', '', ''])
        continue
    n_fwd += 1
    s, e = f['sstart'], f['send']
    hs, he, strand = (s, e, '+') if s <= e else (e, s, '-')
    r = rev.get(blk)
    rbh = 'no'
    rchr = rstart = rend = ''
    if r is not None:
        rs, re_ = r['sstart'], r['send']
        rchr = r['sseqid']
        rstart, rend = (rs, re_) if rs <= re_ else (re_, rs)
        if rchr == qchr and max(qs, rstart) <= min(qe, rend) + MARGIN:
            rbh = 'yes'
            n_rbh += 1
    rows.append([blk, info['trait'], qchr, qs, qe, f['sseqid'], hs, he, strand,
                 f"{f['pident']:.2f}", f"{f['evalue']:.2e}", f"{f['bitscore']:.1f}",
                 rbh, rchr, rstart, rend])

with open(OUT_FILE, 'w', encoding='utf-8') as f:
    f.write('\t'.join(out_header) + '\n')
    for r in rows:
        f.write('\t'.join(str(x) for x in r) + '\n')

print(f"total blocks: {len(blocks)}")
print(f"forward hits: {n_fwd}  ({100.0*n_fwd/len(blocks):.1f}%)")
print(f"RBH confirmed (reciprocal best hit): {n_rbh}  ({100.0*n_rbh/len(blocks):.1f}%)")
print(f"written: {OUT_FILE}")
