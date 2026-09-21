#!/usr/bin/env Rscript
# 01_encode_all_qtls.R
# Block-wise haplotype coding (2 = superior, 0 = inferior) of all candidate QTLs.
# Usage: Rscript 01_encode_all_qtls.R [--base DIR] [--map FILE] [--vcf DIR] [--pheno DIR] [--out DIR]

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

BASE_DIR  <- get_arg("--base",  ".")
MAP_FILE  <- get_arg("--map",   file.path(BASE_DIR, "qtl_block_trait_map.csv"))
VCF_DIR   <- get_arg("--vcf",   file.path(BASE_DIR, "vcf"))
PHENO_DIR <- get_arg("--pheno", file.path(BASE_DIR, "phenotypes"))
OUT_DIR   <- get_arg("--out",   file.path(BASE_DIR, "recode_results"))

# Traits in which higher phenotypic values are agronomically undesirable
REVERSE_TRAITS <- c("VW", "FD")

# Minimum haplotype frequency within a population required to retain a haplotype
MIN_HAP_FREQ <- 0.01

# Population order used for the column order of the coding matrix.
# These labels must match the values in the `bio` column of the phenotype tables.
GROUP_ORDER <- c("Group-A", "Group-B", "Group-C", "Group-D",
                 "Group-E", "Group-F", "Group-G")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("base directory : %s\n", BASE_DIR))
cat(sprintf("QTL-trait map  : %s\n", MAP_FILE))
cat(sprintf("VCF directory  : %s\n", VCF_DIR))
cat(sprintf("phenotype dir  : %s\n", PHENO_DIR))
cat(sprintf("output dir     : %s\n", OUT_DIR))
cat(sprintf("reverse traits : %s\n", paste(REVERSE_TRAITS, collapse = ", ")))
cat(sprintf("min hap freq   : %s\n", MIN_HAP_FREQ))
cat(sprintf("group order    : %s\n", paste(GROUP_ORDER, collapse = ", ")))

# ---- read the QTL-block-trait mapping table ---------------------------------

qtl_map <- read.csv(MAP_FILE, stringsAsFactors = FALSE)
qtl_map$trait_list <- strsplit(qtl_map$Trait_Code, "\\\\")

all_qtls <- sort(unique(qtl_map$QTL_id))
cat(sprintf("QTLs: %d, mapping rows: %d\n", length(all_qtls), nrow(qtl_map)))

# ---- loop over QTLs ---------------------------------------------------------

for (qtl_id in all_qtls) {
  cat("\n========================================\n")
  cat("Processing:", qtl_id, "\n")

  qtl_info   <- qtl_map %>% filter(QTL_id == qtl_id)
  all_traits <- unique(unlist(qtl_info$trait_list))
  cat(sprintf("  traits: %s\n", paste(all_traits, collapse = ", ")))

  qtl_num <- gsub("QTL", "", qtl_id)
  qtl_dir <- file.path(OUT_DIR, paste0("QTL_", qtl_num))
  dir.create(qtl_dir, showWarnings = FALSE, recursive = TRUE)

  # Import each block VCF once and cache the long haplotype table
  all_qtl_blocks <- unique(qtl_info$Block_ID)
  block_hap_cache <- list()
  for (bid in all_qtl_blocks) {
    vcf_file <- file.path(VCF_DIR, paste0(bid, ".vcf"))
    if (!file.exists(vcf_file)) next
    tryCatch({
      vcf <- import_vcf(vcf_file)
      hapResult <- vcf2hap(vcf, hetero_remove = TRUE)
      hap_long <- hapResult[6:nrow(hapResult), ] %>%
        select(Hap, Accession) %>%
        rename(SRR = Accession, haplotype = Hap) %>%
        mutate(block = bid)
      block_hap_cache[[bid]] <- hap_long
    }, error = function(e) {
      cat(sprintf("    warning: failed to import %s: %s\n", bid, e$message))
    })
  }
  cat(sprintf("  loaded %d/%d block haplotype tables\n",
              length(block_hap_cache), length(all_qtl_blocks)))

  # ---- loop over traits -----------------------------------------------------
  for (single_trait in all_traits) {
    cat(sprintf("    trait: %s\n", single_trait))

    is_reverse <- single_trait %in% REVERSE_TRAITS

    # blocks annotated with this trait
    trait_blocks <- qtl_info %>%
      filter(sapply(trait_list, function(tl) single_trait %in% tl)) %>%
      pull(Block_ID)
    trait_blocks <- intersect(trait_blocks, names(block_hap_cache))
    cat(sprintf("      associated blocks: %s\n",
                paste(trait_blocks, collapse = ", ")))
    if (length(trait_blocks) == 0) next

    pheno_file <- file.path(PHENO_DIR, paste0(single_trait, ".csv"))
    if (!file.exists(pheno_file)) {
      cat("      warning: phenotype file not found\n")
      next
    }
    pheno_data <- read.csv(pheno_file, stringsAsFactors = FALSE)
    if (all(is.na(pheno_data[[single_trait]]))) {
      cat("      skipped: all phenotype values are NA\n")
      next
    }

    trait_dir <- file.path(qtl_dir, single_trait)
    dir.create(trait_dir, showWarnings = FALSE, recursive = TRUE)

    qtl_results <- data.frame(
      bio = character(), hap = character(), sio = numeric(),
      sample_count = integer(), stringsAsFactors = FALSE
    )
    qtl_ttest_results <- data.frame(
      bio = character(), block = character(), comparison = character(),
      p_value = numeric(), stringsAsFactors = FALSE
    )
    qtl_sample_haplotypes <- data.frame(
      Accession = character(), block = character(), haplotype = character(),
      stringsAsFactors = FALSE
    )
    qtl_haplotype_sample_counts <- data.frame(
      bio = character(), block = character(), haplotype = character(),
      sample_count = integer(), stringsAsFactors = FALSE
    )

    for (block_id in trait_blocks) {
      hap_long <- block_hap_cache[[block_id]]

      qtl_sample_haplotypes <- rbind(qtl_sample_haplotypes,
        hap_long %>% select(Accession = SRR, block, haplotype))

      merged_data <- hap_long %>% inner_join(pheno_data, by = "SRR")

      for (bio_group in unique(merged_data$bio)) {
        group_data <- merged_data %>% filter(bio == bio_group)
        if (all(is.na(group_data[[single_trait]]))) next

        hap_counts      <- table(group_data$haplotype)
        total_samples   <- nrow(group_data)
        hap_proportions <- hap_counts / total_samples
        valid_haps      <- names(hap_proportions[hap_proportions >= MIN_HAP_FREQ])

        if (length(valid_haps) >= 1) {
          # phenotypic mean and sample size of each retained haplotype
          for (hap_type in valid_haps) {
            hap_pheno <- group_data %>%
              filter(haplotype == hap_type) %>%
              pull(!!single_trait)
            hap_pheno_clean <- hap_pheno[!is.na(hap_pheno)]
            sample_count    <- length(hap_pheno_clean)
            if (sample_count > 0) {
              mean_pheno <- mean(hap_pheno_clean, na.rm = TRUE)
              qtl_results <- rbind(qtl_results, data.frame(
                bio = bio_group,
                hap = paste0(block_id, "-", hap_type),
                sio = mean_pheno,
                sample_count = sample_count,
                stringsAsFactors = FALSE
              ))
              qtl_haplotype_sample_counts <- rbind(
                qtl_haplotype_sample_counts, data.frame(
                  bio = bio_group, block = block_id,
                  haplotype = hap_type, sample_count = sample_count,
                  stringsAsFactors = FALSE
                ))
            }
          }

          # pairwise t-tests between haplotypes of the same block
          if (length(valid_haps) >= 2) {
            hap_pheno_list <- list()
            for (hap_type in valid_haps) {
              hap_data <- group_data %>%
                filter(haplotype == hap_type) %>%
                pull(!!single_trait)
              hap_pheno_list[[hap_type]] <- hap_data[!is.na(hap_data)]
            }
            valid_for_ttest <- sapply(hap_pheno_list, length) >= 2
            if (sum(valid_for_ttest) >= 2) {
              valid_haps_ttest <- valid_haps[valid_for_ttest]
              hap_combinations <- combn(valid_haps_ttest, 2)
              for (i in 1:ncol(hap_combinations)) {
                hap1 <- hap_combinations[1, i]
                hap2 <- hap_combinations[2, i]
                if (length(hap_pheno_list[[hap1]]) >= 2 &&
                    length(hap_pheno_list[[hap2]]) >= 2) {
                  ttest_result <- tryCatch(
                    t.test(hap_pheno_list[[hap1]], hap_pheno_list[[hap2]],
                           alternative = "two.sided"),
                    error = function(e) NULL
                  )
                  if (!is.null(ttest_result)) {
                    qtl_ttest_results <- rbind(qtl_ttest_results, data.frame(
                      bio = bio_group, block = block_id,
                      comparison = paste0(hap1, " vs ", hap2),
                      p_value = ttest_result$p.value,
                      stringsAsFactors = FALSE
                    ))
                  }
                }
              }
            }
          }
        }
      }
    }

    # ---- coding -------------------------------------------------------------
    if (nrow(qtl_results) > 0) {
      if (nrow(qtl_ttest_results) > 0) {
        qtl_ttest_results$significance <- ifelse(
          qtl_ttest_results$p_value < 0.001, "***",
          ifelse(qtl_ttest_results$p_value < 0.01, "**",
                 ifelse(qtl_ttest_results$p_value < 0.05, "*", "ns")))
      }

      qtl_results$block     <- str_extract(qtl_results$hap, "Block_\\d+")
      qtl_results$haplotype <- str_replace(qtl_results$hap, "Block_\\d+-", "")

      qtl_results_encoded <- data.frame()
      for (bg in unique(qtl_results$bio)) {
        for (bn in unique(qtl_results$block)) {
          bd <- qtl_results %>% filter(bio == bg & block == bn)
          if (nrow(bd) > 0) {
            if (nrow(bd) == 1) {
              bd$code <- NA
            } else {
              block_mean <- mean(bd$sio)
              if (is_reverse) {
                # higher values are undesirable: high value -> code 0
                bd$code <- ifelse(bd$sio >= block_mean, 0, 2)
              } else {
                bd$code <- ifelse(bd$sio >= block_mean, 2, 0)
              }
            }
            bd$mean_diff <- bd$sio - mean(bd$sio)
            qtl_results_encoded <- rbind(qtl_results_encoded, bd)
          }
        }
      }
      qtl_results_encoded <- qtl_results_encoded %>%
        select(bio, hap, sio, sample_count, code, mean_diff, block, haplotype)

      # ---- sample x block coding matrix -------------------------------------
      all_samples <- unique(qtl_sample_haplotypes$Accession)
      all_bio_blocks <- expand.grid(
        bio   = unique(qtl_results$bio),
        block = unique(qtl_results$block)
      ) %>% mutate(bio_block = paste(bio, block, sep = "_"))

      all_bio_blocks <- all_bio_blocks %>%
        mutate(
          bio_order = factor(bio, levels = c(GROUP_ORDER,
                                             setdiff(unique(bio), GROUP_ORDER))),
          block_num = as.numeric(str_extract(block, "\\d+"))
        ) %>%
        arrange(bio_order, block_num) %>%
        pull(bio_block)

      sample_code_matrix <- data.frame(Accession = all_samples)

      for (bio_block in all_bio_blocks) {
        bio_name   <- str_split(bio_block, "_Block_")[[1]][1]
        block_name <- paste0("Block_", str_split(bio_block, "_Block_")[[1]][2])

        coding_rules <- qtl_results_encoded %>%
          filter(bio == bio_name & block == block_name) %>%
          select(haplotype, code)

        if (nrow(coding_rules) > 0) {
          block_haps <- qtl_sample_haplotypes %>%
            filter(block == block_name) %>%
            select(Accession, haplotype)
          sample_codes <- block_haps %>%
            left_join(coding_rules, by = "haplotype") %>%
            select(Accession, code) %>%
            rename(!!bio_block := code)
          sample_code_matrix <- sample_code_matrix %>%
            left_join(sample_codes, by = "Accession")
        } else {
          sample_code_matrix[[bio_block]] <- NA
        }
      }

      # ---- write outputs ----------------------------------------------------
      write.csv(qtl_results_encoded,
        file.path(trait_dir, paste0(single_trait, "_haplotype_phenotype_means.csv")),
        row.names = FALSE)
      write.csv(qtl_ttest_results,
        file.path(trait_dir, paste0(single_trait, "_ttest_results.csv")),
        row.names = FALSE)
      write.csv(sample_code_matrix,
        file.path(trait_dir, paste0(single_trait, "_coding_matrix.csv")),
        row.names = FALSE)
      write.csv(qtl_haplotype_sample_counts,
        file.path(trait_dir, paste0(single_trait, "_haplotype_sample_counts.csv")),
        row.names = FALSE)

      if (nrow(qtl_ttest_results) > 0) {
        significance_summary <- qtl_ttest_results %>%
          group_by(bio, block) %>%
          summarise(significant_comparisons = sum(significance != "ns"),
                    total_comparisons = n(), .groups = "drop")
      } else {
        significance_summary <- data.frame(
          bio = character(), block = character(),
          significant_comparisons = integer(),
          total_comparisons = integer(), stringsAsFactors = FALSE)
      }
      write.csv(significance_summary,
        file.path(trait_dir, paste0(single_trait, "_significance_summary.csv")),
        row.names = FALSE)

      cat(sprintf("      done: %d blocks, %d samples\n",
                  length(trait_blocks), nrow(sample_code_matrix)))
    } else {
      cat("      skipped: no valid data\n")
    }
  }
}

cat("\n=== coding complete ===\n")
cat(sprintf("output directory: %s\n", OUT_DIR))
