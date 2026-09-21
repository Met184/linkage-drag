#!/usr/bin/env Rscript
# 05_homeolog_haplotype.R
# Haplotype-phenotype association of the homeologous intervals.
# Usage: Rscript 05_homeolog_haplotype.R [--home DIR] [--ref DIR] [--vcf DIR] [--rbh FILE] [--out DIR]

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

HOME_DIR <- get_arg("--home", ".")
REF_DIR  <- get_arg("--ref",  HOME_DIR)
VCF_DIR  <- get_arg("--vcf",  file.path(HOME_DIR, "block"))
OUT_DIR  <- get_arg("--out",  file.path(HOME_DIR, "haplotype_analysis"))
RBH_FILE <- get_arg("--rbh",  file.path(HOME_DIR, "homeolog_rbh_results.txt"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# Population order used for the row order of the output tables.
BIO_ORDER      <- c("Group-A", "Group-B", "Group-C", "Group-D",
                    "Group-E", "Group-F", "Group-G")
REVERSE_TRAITS <- c("VW", "FD")   # higher values are worse: coding is flipped

# ---------- 1. block -> trait mapping (keep RBH=yes blocks with a VCF) ----------
orth <- read.delim(RBH_FILE, stringsAsFactors = FALSE)
vcf_files <- list.files(VCF_DIR, pattern = "Block_.*_L\\.vcf$", full.names = TRUE)
have_vcf  <- str_extract(basename(vcf_files), "Block_\\d+")
orth <- orth %>% filter(Block_ID %in% have_vcf & RBH == "yes")
cat(sprintf("homeologous blocks to process: %d\n", nrow(orth)))

# ---------- 3. loop over blocks ----------
mean_rows  <- list()
ttest_rows <- list()
count_rows <- list()

for (i in seq_len(nrow(orth))) {
  bid       <- orth$Block_ID[i]
  qchr      <- orth$query_chr[i]
  trait_str <- orth$Trait[i]
  traits    <- trimws(unlist(strsplit(trait_str, "[;/\\\\]")))
  vcf_file  <- file.path(VCF_DIR, paste0(bid, "_L.vcf"))

  cat(sprintf("[%d/%d] %s (%s) trait=%s\n", i, nrow(orth), bid, qchr, trait_str))

  # skip VCFs with 0 SNPs (header only, no variant lines)
  nvar <- sum(!startsWith(readLines(vcf_file, warn = FALSE), "#"))
  if (nvar == 0) { cat("  skipped: 0 SNP\n"); next }

  # import haplotypes
  hap_long <- NULL
  tryCatch({
    vcf <- import_vcf(vcf_file)
    hapResult <- vcf2hap(vcf, hetero_remove = TRUE)
    hap_long <- hapResult[6:nrow(hapResult), ] %>%
      select(Hap, Accession) %>%
      rename(SRR = Accession, haplotype = Hap) %>%
      mutate(block = bid)
  }, error = function(e) {
    cat("  warning: VCF import failed:", conditionMessage(e), "\n")
  })
  if (is.null(hap_long) || nrow(hap_long) == 0) {
    cat("  skipped: no haplotype data\n")
    next
  }

  for (tr in traits) {
    pheno_file <- file.path(REF_DIR, "bio", paste0(tr, ".csv"))
    if (!file.exists(pheno_file)) { cat(sprintf("  skipped: no phenotype file %s\n", tr)); next }
    pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)

    merged <- hap_long %>% inner_join(pheno, by = "SRR")
    if (nrow(merged) == 0 || all(is.na(merged[[tr]]))) next

    for (bg in unique(merged$bio)) {
      gd <- merged %>% filter(bio == bg)
      if (all(is.na(gd[[tr]]))) next

      hap_counts <- table(gd$haplotype)
      valid <- names(hap_counts)[hap_counts / nrow(gd) >= 0.01]
      if (length(valid) < 1) next

      # haplotype phenotype means and coding
      hap_mean <- sapply(valid, function(h) {
        v <- gd %>% filter(haplotype == h) %>% pull(!!tr)
        v <- v[!is.na(v)]
        mean(v)
      })
      hap_n <- sapply(valid, function(h) {
        v <- gd %>% filter(haplotype == h) %>% pull(!!tr)
        sum(!is.na(v))
      })

      block_mean <- mean(hap_mean)
      code <- if (tr %in% REVERSE_TRAITS) {
        ifelse(hap_mean >= block_mean, 0, 2)
      } else {
        ifelse(hap_mean >= block_mean, 2, 0)
      }
      if (length(valid) == 1) code <- NA_integer_

      for (j in seq_along(valid)) {
        h <- valid[j]
        mean_rows[[length(mean_rows) + 1]] <- data.frame(
          block = bid, chr = qchr, trait = tr, bio = bg, haplotype = h,
          mean_pheno = hap_mean[j], sample_count = hap_n[j],
          code = code[j], stringsAsFactors = FALSE)
        count_rows[[length(count_rows) + 1]] <- data.frame(
          block = bid, trait = tr, bio = bg, haplotype = h,
          sample_count = hap_n[j], stringsAsFactors = FALSE)
      }

      # pairwise t-tests (each haplotype needs >= 2 samples)
      if (length(valid) >= 2) {
        hap_list <- lapply(valid, function(h) {
          v <- gd %>% filter(haplotype == h) %>% pull(!!tr)
          v[!is.na(v)]
        })
        names(hap_list) <- valid
        ok <- sapply(hap_list, length) >= 2
        if (sum(ok) >= 2) {
          vh <- valid[ok]
          comb <- combn(vh, 2)
          for (k in seq_len(ncol(comb))) {
            h1 <- comb[1, k]; h2 <- comb[2, k]
            tt <- tryCatch(t.test(hap_list[[h1]], hap_list[[h2]]),
                           error = function(e) NULL)
            if (!is.null(tt)) {
              p <- tt$p.value
              sig <- ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**",
                     ifelse(p < 0.05, "*", "ns")))
              ttest_rows[[length(ttest_rows) + 1]] <- data.frame(
                block = bid, chr = qchr, trait = tr, bio = bg,
                comparison = paste0(h1, " vs ", h2),
                p_value = p, significance = sig, stringsAsFactors = FALSE)
            }
          }
        }
      }
    }
  }
}

# ---------- 4. summarise and write ----------
means  <- if (length(mean_rows) > 0) bind_rows(mean_rows) else data.frame()
ttests <- if (length(ttest_rows) > 0) bind_rows(ttest_rows) else data.frame()
counts <- if (length(count_rows) > 0) bind_rows(count_rows) else data.frame()

if (nrow(means) > 0) {
  means <- means %>% arrange(block, trait, factor(bio, levels = BIO_ORDER), haplotype)
  write.csv(means, file.path(OUT_DIR, "homeolog_haplotype_means.csv"), row.names = FALSE)
}
if (nrow(ttests) > 0) {
  ttests <- ttests %>% arrange(block, trait, factor(bio, levels = BIO_ORDER), p_value)
  write.csv(ttests, file.path(OUT_DIR, "homeolog_ttest_results.csv"), row.names = FALSE)
}
if (nrow(counts) > 0) {
  counts <- counts %>% arrange(block, trait, factor(bio, levels = BIO_ORDER), haplotype)
  write.csv(counts, file.path(OUT_DIR, "homeolog_haplotype_counts.csv"), row.names = FALSE)
}

# block x trait summary: which homeologous intervals show a haplotype-trait association
if (nrow(ttests) > 0) {
  summary <- ttests %>%
    mutate(is_sig = significance != "ns") %>%
    group_by(block, chr, trait) %>%
    summarise(
      n_bio           = n_distinct(bio),
      n_comparisons   = n(),
      n_sig_comparisons = sum(is_sig),
      sig_bio         = paste(sort(unique(bio[is_sig])), collapse = "; "),
      min_p           = min(p_value),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sig_comparisons), min_p, block, trait)
  write.csv(summary, file.path(OUT_DIR, "homeolog_block_summary.csv"), row.names = FALSE)
} else {
  summary <- data.frame()
}

cat("\n=== done ===\n")
cat(sprintf("blocks with haplotype data: %d\n", n_distinct(means$block)))
cat(sprintf("haplotype-phenotype mean rows: %d\n", nrow(means)))
cat(sprintf("t-test rows: %d\n", nrow(ttests)))
cat(sprintf("block x trait with at least one significant comparison: %d\n",
            sum(summary$n_sig_comparisons > 0)))
cat(sprintf("output directory: %s\n", OUT_DIR))
