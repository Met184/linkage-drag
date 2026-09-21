#!/bin/bash
# 07_compensation_by_group.sh
# The same 2x2 table for all accessions, for cultivated and for semi-wild cotton.
# Usage: BASE=dir bash 07_compensation_by_group.sh

set -u

BASE="${BASE:-$PWD}"
BLOCK_DIR="${BLOCK_DIR:-$BASE/block}"
CONSENSUS_DIR="${CONSENSUS_DIR:-$BASE/consensus_coding}"
BIO_DIR="${BIO_DIR:-$BASE/phenotypes}"
OUT_DIR="${OUT_DIR:-$BASE/compensation}"
mkdir -p "$OUT_DIR"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# block|chr|trait|QTL folder|main population
COMBOS=(
"Block_001|Ghir_A01|FM|QTL_001|Group-D"
"Block_002|Ghir_A01|FM|QTL_001|Group-D"
"Block_101|Ghir_A07|FM|QTL_043|Group-A"
"Block_146|Ghir_A07|FS|QTL_054|Group-A"
"Block_155|Ghir_A07|FL|QTL_054|Group-C"
"Block_155|Ghir_A07|FS|QTL_054|Group-C"
"Block_204|Ghir_A08|FS|QTL_071|Group-A"
"Block_269|Ghir_A12|FE|QTL_097|Group-E"
"Block_318|Ghir_D02|LI|QTL_117|Group-B"
"Block_318|Ghir_D02|FWPB|QTL_117|Group-B"
"Block_332|Ghir_D03|FD|QTL_127|Group-B"
"Block_369|Ghir_D04|FE|QTL_139|Group-A"
"Block_372|Ghir_D04|FE|QTL_139|Group-A"
"Block_410|Ghir_D05|FL|QTL_140|Group-A"
"Block_416|Ghir_D05|FL|QTL_142|Group-A"
"Block_570|Ghir_D11|FL|QTL_202|Group-A"
"Block_571|Ghir_D11|FL|QTL_202|Group-A"
"Block_578|Ghir_D11|FL|QTL_203|Group-A"
"Block_605|Ghir_D11|FL|QTL_203|Group-A"
"Block_605|Ghir_D11|FS|QTL_203|Group-A"
)

OUT="$OUT_DIR/compensation_all_pd_results.csv"
{
  echo "block,chr,trait,qtl,main_bio,ALL_n,ALL_SS,ALL_SI,ALL_IS,ALL_II,Cult_n,Cult_SS,Cult_SI,Cult_IS,Cult_II,Semi_n,Semi_SS,Semi_SI,Semi_IS,Semi_II"
} > "$OUT"

for combo in "${COMBOS[@]}"; do
  IFS='|' read -r block chr trait qtl main_bio <<< "$combo"

  vcf="$BLOCK_DIR/${block}_L.vcf"
  consensus="$CONSENSUS_DIR/${qtl}/${trait}_consensus_coding.csv"
  pheno="$BIO_DIR/${trait}.csv"

  [ -f "$vcf" ] || { echo "MISSING VCF: $vcf" >&2; continue; }
  [ -f "$consensus" ] || { echo "MISSING consensus: $consensus" >&2; continue; }
  [ -f "$pheno" ] || { echo "MISSING pheno: $pheno" >&2; continue; }

  reverse=0
  [ "$trait" = "VW" ] || [ "$trait" = "FD" ] && reverse=1

  # 1) VCF -> all samples \t homeolog haplotype (drop samples heterozygous at any SNP)
  awk -F'\t' '
    /^#CHROM/{for(i=10;i<=NF;i++){smp[i]=$i; hap[i]=""; drop[i]=0}; next}
    /^#/{next}
    {
      for(i=10;i<=NF;i++){
        gt=$i; sub(/:.*/,"",gt)
        if(gt=="0/0"||gt=="0|0") hap[i]=hap[i]"0"
        else if(gt=="1/1"||gt=="1|1") hap[i]=hap[i]"1"
        else drop[i]=1
      }
    }
    END{for(i=10;i<=NF;i++) if(!drop[i]) print smp[i]"\t"hap[i]}
  ' "$vcf" > "$TMP/hap.tsv"

  # 2) consensus coding -> all samples \t original code \t pd (0/2 and non-empty pd only)
  awk -F',' -v blk="$block" '
    NR==1{for(i=1;i<=NF;i++){if($i==blk)c=i; if($i=="pd")p=i}; if(c=="")c=-1; next}
    {v=$c; pd=$p; if((v=="0"||v=="2") && pd!="") print $1"\t"v"\t"pd}
  ' "$consensus" > "$TMP/orig.tsv"

  # 3) phenotype -> main-population sample \t trait value
  awk -F',' -v tr="$trait" -v bio="$main_bio" '
    NR==1{for(i=1;i<=NF;i++){if($i=="SRR")s=i; if($i=="bio")b=i; if($i==tr)t=i}; next}
    {if($b==bio && $t!="" && $t!="NA") print $s"\t"$t}
  ' "$pheno" > "$TMP/pheno.tsv"

  # 4) determine the superior/inferior direction of each homeolog haplotype from
  #    the main population; block_mean = mean of haplotype means (unweighted, as
  #    in the R scripts)
  awk -F'\t' -v reverse="$reverse" '
    NR==FNR{ sh[$1]=$2; next }      # hap.tsv: sample -> raw_hap (all samples)
    { if($1 in sh){ h=sh[$1]; v=$2+0; hs[h]+=v; hc[h]++ } }  # pheno.tsv: main population
    END{
      n=0; bm=0
      for(h in hc){ bm += hs[h]/hc[h]; n++ }
      bm = (n>0)? bm/n : 0
      for(h in hc){
        hm = hs[h]/hc[h]
        code = reverse ? (hm<bm?2:0) : (hm>=bm?2:0)
        print h"\t"code
      }
    }
  ' "$TMP/hap.tsv" "$TMP/pheno.tsv" > "$TMP/hap_code.tsv"

  # 5) apply the direction to all samples -> sample \t homeo_code
  awk -F'\t' 'NR==FNR{code[$1]=$2; next} {if($2 in code) print $1"\t"code[$2]}' \
    "$TMP/hap_code.tsv" "$TMP/hap.tsv" > "$TMP/homeo.tsv"

  # 6) merge homeo_code + orig_code + pd
  awk -F'\t' 'NR==FNR{oc[$1]=$2; pd[$1]=$3; next} {if($1 in oc) print $1"\t"$2"\t"oc[$1]"\t"pd[$1]}' \
    "$TMP/orig.tsv" "$TMP/homeo.tsv" > "$TMP/final.tsv"

  # 7) contingency tables: ALL / Cultivated / Semi-wild
  awk -F'\t' -v block="$block" -v chr="$chr" -v trait="$trait" -v qtl="$qtl" -v main_bio="$main_bio" '
    {
      h=$2; o=$3; p=$4
      if(o==2 && h==2){ SS++; if(p=="Cultivated")cSS++; if(p=="Semi-wild")sSS++ }
      else if(o==2 && h==0){ SI++; if(p=="Cultivated")cSI++; if(p=="Semi-wild")sSI++ }
      else if(o==0 && h==2){ IS++; if(p=="Cultivated")cIS++; if(p=="Semi-wild")sIS++ }
      else if(o==0 && h==0){ II++; if(p=="Cultivated")cII++; if(p=="Semi-wild")sII++ }
    }
    END{
      printf "%s,%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
        block, chr, trait, qtl, main_bio,
        SS+SI+IS+II, SS, SI, IS, II,
        cSS+cSI+cIS+cII, cSS, cSI, cIS, cII,
        sSS+sSI+sIS+sII, sSS, sSI, sIS, sII
    }
  ' "$TMP/final.tsv" >> "$OUT"
done

echo "DONE -> $OUT"
cat "$OUT"
