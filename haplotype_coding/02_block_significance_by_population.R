#!/usr/bin/env Rscript
# 02_block_significance_by_population.R
# Significance and effect direction per QTL x block x trait and population.
# Usage: Rscript 02_block_significance_by_population.R [--coding DIR] [--out DIR]

suppressMessages({
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

WORK     <- get_arg("--work", ".")
CODE_DIR <- get_arg("--coding", file.path(WORK, "coding"))
OUT_DIR  <- get_arg("--out", file.path(WORK, "block_qc"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ALL_GROUPS <- c("Group-A", "Group-B", "Group-C", "Group-D",
                "Group-E", "Group-F", "Group-G")

# ---- scan all QTL / trait directories ----
qtl_dirs <- list.dirs(CODE_DIR, recursive = FALSE, full.names = TRUE)
cat(sprintf("found %d QTL directories\n", length(qtl_dirs)))

all_rows <- list()

for (qtl_dir in qtl_dirs) {
  qtl_name <- basename(qtl_dir)
  trait_dirs <- list.dirs(qtl_dir, recursive = FALSE, full.names = TRUE)

  for (trait_dir in trait_dirs) {
    trait_name <- basename(trait_dir)

    means_file <- file.path(trait_dir, paste0(trait_name, "_haplotype_phenotype_means.csv"))
    ttest_file <- file.path(trait_dir, paste0(trait_name, "_ttest_results.csv"))

    if (!file.exists(means_file)) next

    means <- read.csv(means_file, stringsAsFactors = FALSE)

    if (file.exists(ttest_file)) {
      ttest <- read.csv(ttest_file, stringsAsFactors = FALSE)
    } else {
      ttest <- data.frame(bio = character(), block = character(),
                          comparison = character(), p_value = numeric(),
                          significance = character(), stringsAsFactors = FALSE)
    }

    blocks <- unique(means$block)

    for (blk in blocks) {
      blk_means <- means %>% filter(block == blk)

      row_data <- list(
        QTL   = qtl_name,
        Block = blk,
        Trait = trait_name
      )

      hap2_sets <- list()      # code-2 haplotype sets per population
      hap_code_map <- list()   # haplotype -> list(population = code)

      for (grp in ALL_GROUPS) {
        grp_data <- blk_means %>% filter(bio == grp)

        if (nrow(grp_data) == 0) {
          row_data[[grp]] <- "-"
          next
        }

        hap2 <- sort(grp_data$haplotype[grp_data$code == 2])
        hap0 <- sort(grp_data$haplotype[grp_data$code == 0])

        hap2_str <- paste(hap2, collapse = ",")
        hap2_sets[[grp]] <- hap2_str

        # record the code of each haplotype in this population
        for (j in 1:nrow(grp_data)) {
          h <- grp_data$haplotype[j]
          c <- grp_data$code[j]
          if (is.null(hap_code_map[[h]])) hap_code_map[[h]] <- list()
          hap_code_map[[h]][[grp]] <- ifelse(is.na(c), NA, c)
        }

        # direction of the difference between code 2 and code 0
        if (length(hap2) > 0 && length(hap0) > 0) {
          mean2 <- mean(grp_data$sio[grp_data$code == 2])
          mean0 <- mean(grp_data$sio[grp_data$code == 0])
          diff_val <- mean2 - mean0
          dir_sign <- ifelse(diff_val > 0, "+", "-")
        } else {
          dir_sign <- "?"
        }

        # significance: best (lowest) level among the comparisons of this block
        grp_ttest <- ttest %>% filter(bio == grp & block == blk)
        if (nrow(grp_ttest) > 0 && any(grp_ttest$significance %in% c("***", "**", "*"))) {
          sig_levels <- grp_ttest$significance
          if ("***" %in% sig_levels) sig <- "***"
          else if ("**" %in% sig_levels) sig <- "**"
          else sig <- "*"
        } else if (nrow(grp_ttest) > 0) {
          sig <- "ns"
        } else {
          sig <- "ns"
        }

        # cell format: significance | code-2 haplotypes | direction
        row_data[[grp]] <- sprintf("%s|%s|%s", sig, hap2_str, dir_sign)
      }

      # ---- direction consistency, judged haplotype by haplotype ----
      consistent_high <- c()  # code 2 in 2 or more populations
      consistent_low  <- c()  # code 0 in 2 or more populations
      mixed_haps      <- c()  # conflicting code between populations

      for (hap in names(hap_code_map)) {
        codes <- unlist(hap_code_map[[hap]])
        codes_valid <- codes[!is.na(codes)]
        if (length(codes_valid) >= 2) {  # present in at least 2 populations
          if (all(codes_valid == 2)) {
            consistent_high <- c(consistent_high, hap)
          } else if (all(codes_valid == 0)) {
            consistent_low <- c(consistent_low, hap)
          } else {
            mixed_haps <- c(mixed_haps, hap)
          }
        }
      }

      row_data$consistent_superior_haplotypes  <- paste(consistent_high, collapse = ",")
      row_data$consistent_inferior_haplotypes  <- paste(consistent_low, collapse = ",")
      row_data$mixed_direction_haplotypes      <- paste(mixed_haps, collapse = ",")

      # a consistently superior or inferior haplotype means the direction agrees
      if (length(consistent_high) + length(consistent_low) > 0) {
        row_data$direction_consistency <- "consistent"
        valid_sets <- hap2_sets[hap2_sets != ""]
        if (length(unique(valid_sets)) > 1) {
          freq <- table(as.character(valid_sets))
          consensus_set <- names(which.max(freq))
          row_data$direction_consistent_groups   <- paste(names(valid_sets)[valid_sets == consensus_set], collapse = ", ")
          row_data$direction_inconsistent_groups <- paste(names(valid_sets)[valid_sets != consensus_set], collapse = ", ")
        } else {
          row_data$direction_consistent_groups   <- paste(names(valid_sets), collapse = ", ")
          row_data$direction_inconsistent_groups <- ""
        }
      } else if (length(mixed_haps) > 0) {
        row_data$direction_consistency <- "inconsistent"
        trouble_grps <- c()
        for (hap in mixed_haps) {
          hap_grp_codes <- hap_code_map[[hap]]
          vals <- unlist(hap_grp_codes)
          grps_2 <- names(vals)[vals == 2]
          grps_0 <- names(vals)[vals == 0]
          trouble_grps <- c(trouble_grps, grps_2, grps_0)
        }
        all_grps_with_data <- names(hap2_sets)[hap2_sets != ""]
        trouble_grps <- unique(trouble_grps)
        row_data$direction_consistent_groups   <- paste(setdiff(all_grps_with_data, trouble_grps), collapse = ", ")
        row_data$direction_inconsistent_groups <- paste(trouble_grps, collapse = ", ")
      } else {
        row_data$direction_consistency <- ifelse(length(hap2_sets[hap2_sets != ""]) >= 2, "single_group", "no_data")
        row_data$direction_consistent_groups   <- paste(names(hap2_sets)[hap2_sets != ""], collapse = ", ")
        row_data$direction_inconsistent_groups <- ""
      }

      # ---- counts ----
      sig_count <- 0
      for (grp in ALL_GROUPS) {
        val <- row_data[[grp]]
        if (val != "-") {
          sig_part <- strsplit(val, "\\|")[[1]][1]
          if (sig_part %in% c("***", "**", "*")) sig_count <- sig_count + 1
        }
      }
      row_data$n_significant_groups <- sig_count
      row_data$n_groups_with_data   <- sum(sapply(ALL_GROUPS, function(g) row_data[[g]] != "-"))

      # ---- per-population code-2 haplotype list (plain text, for filtering) ----
      for (grp in ALL_GROUPS) {
        col_name <- paste0(grp, "_hap2")
        row_data[[col_name]] <- if (grp %in% names(hap2_sets)) hap2_sets[[grp]] else ""
      }

      # ---- per-population mean(code 2) - mean(code 0) ----
      for (grp in ALL_GROUPS) {
        col_name <- paste0(grp, "_diff")
        grp_data <- blk_means %>% filter(bio == grp)
        if (nrow(grp_data) > 0) {
          m2 <- mean(grp_data$sio[grp_data$code == 2])
          m0 <- mean(grp_data$sio[grp_data$code == 0])
          if (!is.na(m2) && !is.na(m0)) {
            row_data[[col_name]] <- round(m2 - m0, 4)
          } else {
            row_data[[col_name]] <- NA
          }
        } else {
          row_data[[col_name]] <- NA
        }
      }

      all_rows[[length(all_rows) + 1]] <- as.data.frame(row_data, stringsAsFactors = FALSE)
    }
  }
}

# ---- merge and write ----
if (length(all_rows) > 0) {
  result <- bind_rows(all_rows)

  aux_cols <- c(
    paste0(ALL_GROUPS, "_hap2"),
    paste0(ALL_GROUPS, "_diff")
  )
  main_cols <- c("QTL", "Block", "Trait", ALL_GROUPS,
                 "direction_consistency", "direction_consistent_groups",
                 "direction_inconsistent_groups",
                 "consistent_superior_haplotypes", "consistent_inferior_haplotypes",
                 "mixed_direction_haplotypes",
                 "n_significant_groups", "n_groups_with_data")
  all_cols <- c(main_cols, aux_cols)
  all_cols <- intersect(all_cols, names(result))
  result <- result[, all_cols]

  write.csv(result, file.path(OUT_DIR, "block_significance_summary.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  cat(sprintf("\n=== done ===\n"))
  cat(sprintf("rows: %d\n", nrow(result)))
  cat(sprintf("QTLs: %d\n", n_distinct(result$QTL)))
  cat(sprintf("blocks: %d\n", n_distinct(result$Block)))
  cat(sprintf("traits: %d\n", n_distinct(result$Trait)))
  cat(sprintf("\ndirection consistency:\n"))
  print(table(result$direction_consistency))
  cat(sprintf("\nnumber of significant populations:\n"))
  print(table(result$n_significant_groups))
  cat(sprintf("\nwritten to: %s\n", file.path(OUT_DIR, "block_significance_summary.csv")))
} else {
  cat("warning: no data found\n")
}
