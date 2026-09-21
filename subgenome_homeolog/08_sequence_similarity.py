#!/usr/bin/env python3
"""08_sequence_similarity.py

Sequence similarity of the homeologous intervals, block level and per sample.

Usage: python3 08_sequence_similarity.py --work DIR --genome GENOME.fasta --master-vcf FILE
"""

import argparse
import os
import csv
import subprocess
from collections import Counter, defaultdict
import difflib

# ============ configuration ============
ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--work', default='.', help='working directory (default: current directory)')
ap.add_argument('--genome', default='genome.fasta',
                help='reference genome FASTA, indexed with samtools faidx')
ap.add_argument('--master-vcf', default='all_samples.vcf.gz',
                help='VCF containing all samples, used to extract homeologous sequences')
ap.add_argument('--res-dir', default=None,
                help='directory with the block haplotype results and block VCFs '
                     '(default: <work>/491block)')
ap.add_argument('--out', default=None,
                help='output directory (default: <work>/sequence_similarity)')
args = ap.parse_args()

WORK       = args.work
GENOME     = args.genome
MASTER_VCF = args.master_vcf
QTL_FILE   = os.path.join(WORK, "QTL_with_topSNP.txt")
RBH_FILE   = os.path.join(WORK, "homeolog_rbh_results.txt")
RES_DIR    = args.res_dir if args.res_dir else os.path.join(WORK, "491block")
OUT_DIR    = args.out if args.out else os.path.join(WORK, "sequence_similarity")

MEANS_CSV       = os.path.join(RES_DIR, "haplotype_means_coding.csv")
TTEST_CSV       = os.path.join(RES_DIR, "haplotype_ttest_results.csv")
SAMPLE_CODE_CSV = os.path.join(RES_DIR, "haplotype_sample_code.csv")

os.makedirs(OUT_DIR, exist_ok=True)

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def revcomp(s):
    return s.translate(COMP)[::-1]

# ============ basic IO ============
def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def read_csv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))

def faidx(region):
    """Fetch a sequence with samtools faidx; return the uppercase sequence or None."""
    r = subprocess.run(["samtools", "faidx", GENOME, region],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None
    return "".join(l.strip() for l in r.stdout.splitlines() if not l.startswith(">")).upper()

def parse_vcf_lines(lines):
    """Parse VCF text lines -> (samples, snps). snps: [{pos,ref,alt,gt:{sample:GT}}]"""
    samples, snps = [], []
    for line in lines:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.strip().split("\t")[9:]
            continue
        if line.startswith("#"):
            continue
        c = line.strip().split("\t")
        if len(c) < 10:
            continue
        pos = int(c[1]); ref = c[3].upper(); alt = c[4].split(",")[0].upper()
        gt = {}
        for i, s in enumerate(samples):
            gt[s] = c[9 + i].split(":")[0]
        snps.append({"pos": pos, "ref": ref, "alt": alt, "gt": gt})
    return samples, snps

def parse_vcf(vcf_file):
    if not os.path.exists(vcf_file):
        return None, None
    with open(vcf_file, encoding="utf-8") as f:
        return parse_vcf_lines(f)

def vcf_region_lines(region):
    """Extract a region from the master VCF with `bcftools view -r`; None on failure."""
    r = subprocess.run(["bcftools", "view", "-r", region, MASTER_VCF],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None
    return r.stdout.splitlines()

def map_positions(orig_ref, homeo_ref):
    """Globally align the two references and return an original 1-based ->
    homeolog 1-based map built from equal/replace blocks."""
    sm = difflib.SequenceMatcher(None, orig_ref, homeo_ref, autojunk=False)
    mapping = {}
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag in ("equal", "replace"):
            n = min(i2 - i1, j2 - j1)
            for k in range(n):
                mapping[i1 + k + 1] = j1 + k + 1
    return mapping

def consensus_allele(gts):
    """Majority allele among a set of sample genotypes: 'REF' / 'ALT' / None."""
    c = Counter(gts)
    n_ref = c.get("0/0", 0) + c.get("0|0", 0)
    n_alt = c.get("1/1", 0) + c.get("1|1", 0)
    if n_ref == 0 and n_alt == 0:
        return None
    return "REF" if n_ref >= n_alt else "ALT"

# ============ 1. read QTL / RBH tables ============
qtl = {}
for r in read_tsv(QTL_FILE):
    bid = r["Block_ID"].strip()
    if not bid or not r["new_start"].strip():
        continue
    qtl[bid] = {"chr": r["chr"].strip(), "trait": r["Trait"].strip(),
                "start": int(r["new_start"]), "end": int(r["new_end"])}

rbh = {}
for r in read_tsv(RBH_FILE):
    bid = r["Block_ID"].strip()
    hc = r["homeolog_chr"].strip()
    rbh[bid] = {
        "homeolog_chr": hc,
        "homeolog_start": int(r["homeolog_start"]) if hc else None,
        "homeolog_end": int(r["homeolog_end"]) if hc else None,
        "strand": r["strand"].strip(),
        "RBH": r["RBH"].strip(),
    }

# ============ 2. keep significant and direction-consistent blocks ============
sig_blocks = defaultdict(int)
for r in read_csv(TTEST_CSV):
    if r["significance"] in ("*", "**", "***"):
        sig_blocks[r["block"]] += 1

hap_codes = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
for r in read_csv(MEANS_CSV):
    if r["code"] == "2":
        hap_codes[r["block"]][r["haplotype"]]["sup"] += 1
    elif r["code"] == "0":
        hap_codes[r["block"]][r["haplotype"]]["inf"] += 1

def direction_consistent(bid):
    for h, d in hap_codes[bid].items():
        if d.get("sup", 0) > 0 and d.get("inf", 0) > 0:
            return False
    return True

sample_code = defaultdict(dict)
for r in read_csv(SAMPLE_CODE_CSV):
    sample_code[r["block"]][r["SRR"]] = r["code"]

valid_blocks = [bid for bid in qtl
                if bid in sig_blocks and direction_consistent(bid)]
valid_blocks.sort()

print(f"total QTL blocks: {len(qtl)}")
print(f"significant blocks (>=1 significant t-test): {len(sig_blocks)}")
print(f"significant and direction-consistent blocks: {len(valid_blocks)}")

# ============ 3. loop over blocks ============
rows = []
per_sample_rows = []
fa_sup_inf = []
fa_homeo = []
log = []

for bid in valid_blocks:
    q = qtl[bid]
    qchr, trait = q["chr"], q["trait"]
    rb = rbh.get(bid, {})
    hchr = rb.get("homeolog_chr", "")
    strand = rb.get("strand", "+")
    is_rbh = (rb.get("RBH", "") == "yes")

    try:
        orig_samples, orig_snps = parse_vcf(os.path.join(RES_DIR, bid + ".vcf"))
        if not orig_snps:
            log.append(f"{bid}: no SNP in the original VCF, skipped"); continue
        sup_srrs = [s for s, c in sample_code.get(bid, {}).items() if c == "2"]
        inf_srrs = [s for s, c in sample_code.get(bid, {}).items() if c == "0"]
        if not sup_srrs or not inf_srrs:
            log.append(f"{bid}: missing superior or inferior samples, skipped"); continue
        sup_set = set(sup_srrs); inf_set = set(inf_srrs)

        orig_ref = faidx(f"{qchr}:{q['start']}-{q['end']}")
        if orig_ref is None or len(orig_ref) != (q["end"] - q["start"] + 1):
            log.append(f"{bid}: failed to extract the original 201 bp sequence, skipped"); continue

        # For every significant SNP take the majority allele in the superior and
        # inferior groups; keep only SNPs that discriminate the two groups.
        disc = []  # (orig_rel_1based, sup_base, inf_base)
        seq_sup = list(orig_ref); seq_inf = list(orig_ref)
        for snp in orig_snps:
            orig_rel = snp["pos"] - q["start"] + 1
            if orig_rel < 1 or orig_rel > len(orig_ref):
                continue
            a_sup = consensus_allele(snp["gt"][s] for s in sup_set if s in snp["gt"])
            a_inf = consensus_allele(snp["gt"][s] for s in inf_set if s in snp["gt"])
            if a_sup is None or a_inf is None or a_sup == a_inf:
                continue
            b_sup = snp["ref"] if a_sup == "REF" else snp["alt"]
            b_inf = snp["ref"] if a_inf == "REF" else snp["alt"]
            disc.append((orig_rel, b_sup, b_inf))
            seq_sup[orig_rel - 1] = b_sup
            seq_inf[orig_rel - 1] = b_inf

        if not disc:
            log.append(f"{bid}: no SNP discriminating superior from inferior, skipped"); continue

        fa_sup_inf.append(f">{bid}|superior|{trait}\n{''.join(seq_sup)}\n")
        fa_sup_inf.append(f">{bid}|inferior|{trait}\n{''.join(seq_inf)}\n")

        row = {
            "block": bid, "chr": qchr, "trait": trait,
            "n_sig_comparisons": sig_blocks[bid],
            "n_snps": len(orig_snps),
            "n_discriminating": len(disc),
            "n_sup_samples": len(sup_srrs),
            "n_inf_samples": len(inf_srrs),
            "homeolog_chr": hchr,
            "homeolog_start": rb.get("homeolog_start", ""),
            "homeolog_end": rb.get("homeolog_end", ""),
            "strand": strand,
            "RBH": "yes" if is_rbh else "no",
        }

        # --- homeologous interval: block-level reference + per-sample master VCF ---
        if hchr:
            homeo_ref = faidx(f"{hchr}:{rb['homeolog_start']}-{rb['homeolog_end']}")
            if homeo_ref is None:
                log.append(f"{bid}: failed to extract the homeologous sequence")
            else:
                aligned = revcomp(homeo_ref) if strand == "-" else homeo_ref
                mapping = map_positions(orig_ref, aligned)

                # block level: homeologous reference sequence vs superior/inferior
                score = match_sup = match_inf = n_mapped = 0
                for orig_rel, bs, bi in disc:
                    hr = mapping.get(orig_rel)
                    if hr is None:
                        continue
                    hb = aligned[hr - 1]
                    n_mapped += 1
                    if hb == bs:
                        score += 1; match_sup += 1
                    elif hb == bi:
                        score -= 1; match_inf += 1
                homeo_cls = ("superior-like" if score > 0 else
                             ("inferior-like" if score < 0 else "ambiguous"))
                row.update({
                    "homeo_class": homeo_cls, "score": score,
                    "n_mapped": n_mapped, "match_sup": match_sup, "match_inf": match_inf,
                })

                # per sample: extract the homeologous region from the master VCF,
                # rebuild each sample's 201 bp sequence and classify it
                region = f"{hchr}:{rb['homeolog_start']}-{rb['homeolog_end']}"
                hlines = vcf_region_lines(region)
                hsamples, hsnps = (parse_vcf_lines(hlines) if hlines else (None, None))
                if hsamples:
                    hstart = rb['homeolog_start']
                    for s in hsamples:
                        seq = list(homeo_ref)
                        for snp in hsnps:
                            rel = snp["pos"] - hstart + 1
                            if rel < 1 or rel > len(homeo_ref):
                                continue
                            gt = snp["gt"].get(s)
                            if gt in ("0/0", "0|0"):
                                seq[rel - 1] = snp["ref"]
                            elif gt in ("1/1", "1|1"):
                                seq[rel - 1] = snp["alt"]
                        aseq = revcomp("".join(seq)) if strand == "-" else "".join(seq)
                        fa_homeo.append(f">{bid}|{s}|{trait}\n{aseq}\n")

                        sc = msup = minf = 0
                        for orig_rel, bs, bi in disc:
                            hr = mapping.get(orig_rel)
                            if hr is None:
                                continue
                            hb = aseq[hr - 1]
                            if hb == bs:
                                sc += 1; msup += 1
                            elif hb == bi:
                                sc -= 1; minf += 1
                        cls = ("superior-like" if sc > 0 else
                               ("inferior-like" if sc < 0 else "ambiguous"))
                        per_sample_rows.append({
                            "block": bid, "chr": qchr, "trait": trait, "SRR": s,
                            "orig_code": sample_code.get(bid, {}).get(s, ""),
                            "homeo_class": cls, "score": sc,
                            "match_sup": msup, "match_inf": minf,
                        })
        else:
            row["homeo_class"] = "no_homeolog"

        rows.append(row)
    except Exception as e:
        log.append(f"{bid}: exception {e}")

# ============ 4. write results ============
if rows:
    cols = ["block", "chr", "trait", "n_sig_comparisons", "n_snps",
            "n_discriminating", "n_sup_samples", "n_inf_samples",
            "homeolog_chr", "homeolog_start", "homeolog_end", "strand", "RBH",
            "homeo_class", "score", "n_mapped", "match_sup", "match_inf"]
    with open(os.path.join(OUT_DIR, "homeolog_classification.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

if per_sample_rows:
    pcols = ["block", "chr", "trait", "SRR", "orig_code",
             "homeo_class", "score", "match_sup", "match_inf"]
    with open(os.path.join(OUT_DIR, "homeolog_per_sample_classification.csv"),
              "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=pcols, extrasaction="ignore")
        w.writeheader(); w.writerows(per_sample_rows)

with open(os.path.join(OUT_DIR, "superior_inferior_sequences.fa"), "w", encoding="utf-8") as f:
    f.write("".join(fa_sup_inf))

with open(os.path.join(OUT_DIR, "homeolog_sample_sequences.fa"), "w", encoding="utf-8") as f:
    f.write("".join(fa_homeo))

with open(os.path.join(OUT_DIR, "run_summary.txt"), "w", encoding="utf-8") as f:
    f.write(f"total QTL blocks: {len(qtl)}\n")
    f.write(f"significant blocks: {len(sig_blocks)}\n")
    f.write(f"significant and direction-consistent blocks: {len(valid_blocks)}\n")
    f.write(f"blocks with discriminating SNPs and output: {len(rows)}\n")
    if rows:
        cls_cnt = Counter(r.get("homeo_class", "") for r in rows)
        for k, v in cls_cnt.items():
            f.write(f"  {k}: {v}\n")
        n_rbh = sum(1 for r in rows if r["RBH"] == "yes")
        n_fwd_only = sum(1 for r in rows if r["RBH"] != "yes" and r.get("homeolog_chr", ""))
        n_no_homolog = sum(1 for r in rows if not r.get("homeolog_chr", ""))
        f.write(f"RBH=yes: {n_rbh}  forward only: {n_fwd_only}  no_homeolog: {n_no_homolog}\n")
    f.write(f"per-sample homeologous sequences: {len(fa_homeo)}\n")
    f.write(f"per-sample classification rows: {len(per_sample_rows)}\n")
    f.write("\n--- skipped / exception log ---\n")
    f.write("\n".join(log) + "\n")

print(f"blocks with discriminating SNPs and output: {len(rows)}")
print(f"per-sample classification rows: {len(per_sample_rows)}")
print(f"output directory: {OUT_DIR}")
