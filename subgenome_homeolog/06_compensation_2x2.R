#!/usr/bin/env Rscript
# 06_compensation_2x2.R
# 2x2 contingency table between an original block and its homeologous interval.
# Usage: Rscript 06_compensation_2x2.R [--home DIR] [--consensus DIR] [--ref DIR] [--out DIR]

suppressMessages({
  library(geneHapR)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ---- CONFIG -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1L) args[i + 1L] else default
}

HOME_DIR      <- get_arg("--home", ".")
CONSENSUS_DIR <- get_arg("--consensus", file.path(HOME_DIR, "consensus_coding"))
REF_DIR       <- get_arg("--ref", HOME_DIR)
VCF_DIR       <- get_arg("--vcf", file.path(HOME_DIR, "block"))
OUT_DIR       <- get_arg("--out", file.path(HOME_DIR, "compensation"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


REVERSE_TRAITS <- c("VW", "FD")

# ---------- block x trait combinations analysed here (block, chr, trait, QTL folder, main population) ----------
combos <- data.frame(
  block = c("Block_001","Block_002","Block_101","Block_146","Block_155","Block_155",
            "Block_204","Block_269","Block_318","Block_318","Block_332","Block_369",
            "Block_372","Block_410","Block_416","Block_570","Block_571","Block_578",
            "Block_605","Block_605"),
  chr   = c("Ghir_A01","Ghir_A01","Ghir_A07","Ghir_A07","Ghir_A07","Ghir_A07",
            "Ghir_A08","Ghir_A12","Ghir_D02","Ghir_D02","Ghir_D03","Ghir_D04",
            "Ghir_D04","Ghir_D05","Ghir_D05","Ghir_D11","Ghir_D11","Ghir_D11",
            "Ghir_D11","Ghir_D11"),
  trait = c("FM","FM","FM","FS","FL","FS","FS","FE","LI","FWPB","FD","FE",
            "FE","FL","FL","FL","FL","FL","FL","FS"),
  qtl   = c("QTL_001","QTL_001","QTL_043","QTL_054","QTL_054","QTL_054",
            "QTL_071","QTL_097","QTL_117","QTL_117","QTL_127","QTL_139",
            "QTL_139","QTL_140","QTL_142","QTL_202","QTL_202","QTL_203",
            "QTL_203","QTL_203"),
  main_bio = c("Group-D","Group-D","Group-A","Group-A","Group-C","Group-C","Group-A","Group-E",
               "Group-B","Group-B","Group-B","Group-A","Group-A","Group-A","Group-A","Group-A",
               "Group-A","Group-A","Group-A","Group-A"),
  stringsAsFactors = FALSE
)

results <- list()

for (i in seq_len(nrow(combos))) {
  block    <- combos$block[i]
  chr      <- combos$chr[i]
  trait    <- combos$trait[i]
  qtl      <- combos$qtl[i]
  main_bio <- combos$main_bio[i]

  vcf_file       <- file.path(VCF_DIR, paste0(block, "_L.vcf"))
  consensus_file <- file.path(CONSENSUS_DIR, qtl, paste0(trait, "_consensus_coding.csv"))
  pheno_file     <- file.path(REF_DIR, "bio", paste0(trait, ".csv"))

  cat(sprintf("[%d/%d] %s %s %s\n", i, nrow(combos), block, trait, main_bio))
  if (!file.exists(vcf_file) || !file.exists(consensus_file) || !file.exists(pheno_file)) {
    cat("  skipped: missing file\n"); next
  }

  # ---- 1. homeolog sample haplotypes (same as 05_homeolog_haplotype.R) ----
  hap_long <- NULL
  tryCatch({
    vcf <- import_vcf(vcf_file)
    hapResult <- vcf2hap(vcf, hetero_remove = TRUE)
    hap_long <- hapResult[6:nrow(hapResult), ] %>%
      select(Hap, Accession) %>%
      rename(SRR = Accession, haplotype = Hap) %>%
      filter(!is.na(haplotype) & haplotype != "" & haplotype != "NA")
  }, error = function(e) {
    cat("  warning: VCF import failed:", conditionMessage(e), "\n")
  })
  if (is.null(hap_long) || nrow(hap_long) == 0) { cat("  skipped: no haplotype\n"); next }

  # ---- 2. consensus coding of the original block (Accession -> orig_code 0/2) ----
  cons <- read.csv(consensus_file, stringsAsFactors = FALSE)
  if (!block %in% names(cons)) { cat("  skipped: consensus coding lacks this block column\n"); next }
  orig <- cons[, c("Accession", block)]
  names(orig) <- c("SRR", "orig_code")
  orig <- orig[!is.na(orig$orig_code) & orig$orig_code != "" & orig$orig_code != "NA", ]
  orig$orig_code <- as.numeric(orig$orig_code)
  orig <- orig[orig$orig_code %in% c(0, 2), ]

  # ---- 3. phenotype (main population only) ----
  pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
  pheno <- pheno[pheno$bio == main_bio & !is.na(pheno[[trait]]), c("SRR", trait)]
  names(pheno) <- c("SRR", "pheno")

  # ---- 4. merge homeolog haplotype + phenotype -> haplotype coding ----
  hp <- merge(hap_long, pheno, by = "SRR")
  if (nrow(hp) == 0) { cat("  skipped: no phenotype match\n"); next }

  hap_mean_df <- hp %>%
    group_by(haplotype) %>%
    summarise(hap_mean = mean(pheno), .groups = "drop")
  block_mean <- mean(hap_mean_df$hap_mean)   # same as 05_homeolog_haplotype.R (mean of haplotype means)
  if (trait %in% REVERSE_TRAITS) {
    hap_mean_df$code <- ifelse(hap_mean_df$hap_mean >= block_mean, 0, 2)
  } else {
    hap_mean_df$code <- ifelse(hap_mean_df$hap_mean >= block_mean, 2, 0)
  }
  hp <- merge(hp, hap_mean_df[, c("haplotype", "code")], by = "haplotype")

  # ---- 5. merge original coding -> contingency table ----
  df <- merge(hp, orig, by = "SRR")
  if (nrow(df) == 0) { cat("  skipped: no original coding match\n"); next }

  SS <- sum(df$orig_code == 2 & df$code == 2)
  SI <- sum(df$orig_code == 2 & df$code == 0)
  IS <- sum(df$orig_code == 0 & df$code == 2)
  II <- sum(df$orig_code == 0 & df$code == 0)

  n2 <- SS + SI
  n0 <- IS + II
  pct2 <- if (n2 > 0) SS / n2 else NA
  pct0 <- if (n0 > 0) IS / n0 else NA

  # Fisher's exact test
  mat <- matrix(c(SS, SI, IS, II), nrow = 2, byrow = TRUE)
  ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
  fisher_p <- if (!is.null(ft)) ft$p.value else NA
  OR        <- if (!is.null(ft)) as.numeric(ft$estimate) else NA

  rel <- if (is.na(pct2)) "no_superior_original_samples"
         else if (pct2 > 0.5) "concordant_homeolog_superior"
         else if (pct2 < 0.5) "compensated_homeolog_inferior"
         else "no_preference_0.5"

  results[[length(results) + 1]] <- data.frame(
    block = block, chr = chr, trait = trait, qtl = qtl, main_bio = main_bio,
    n_SS = SS, n_SI = SI, n_IS = IS, n_II = II,
    pct_homeo2_given_orig2 = round(pct2, 3),
    pct_homeo2_given_orig0 = round(pct0, 3),
    odds_ratio = round(OR, 3),
    fisher_p = signif(fisher_p, 3),
    relationship = rel,
    stringsAsFactors = FALSE
  )
}

res <- if (length(results) > 0) bind_rows(results) else data.frame()
if (nrow(res) > 0) {
  write.csv(res, file.path(OUT_DIR, "compensation_results_R.csv"), row.names = FALSE)
  cat("\n=== results ===\n")
  print(res, row.names = FALSE)
  cat(sprintf("\nconcordant (homeolog superior): %d\n", sum(res$relationship == "concordant_homeolog_superior")))
  cat(sprintf("compensated (homeolog inferior): %d\n", sum(res$relationship == "compensated_homeolog_inferior")))
  cat(sprintf("no superior original samples: %d\n", sum(res$relationship == "no_superior_original_samples")))
  cat(sprintf("output: %s/compensation_results_R.csv\n", OUT_DIR))
} else {
  cat("no results\n")
}
