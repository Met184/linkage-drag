#!/usr/bin/env Rscript
# cross_qtl_epistasis_analysis.R
# Epistasis between the A07 drag QTLs (QTL051-QTL054) and genome-wide FL/FS QTLs.
# Usage: Rscript cross_qtl_epistasis_analysis.R [--work DIR] [--consensus DIR] [--out DIR]

suppressMessages({
  library(data.table)
  library(dplyr)
})

# ---- CONFIG -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1L) args[i + 1L] else default
}

WORK_DIR      <- get_arg("--work", ".")
CONSENSUS_DIR <- get_arg("--consensus", file.path(WORK_DIR, "consensus_coding"))
OUT_DIR       <- get_arg("--out", WORK_DIR)
BLOCK_FILE    <- get_arg("--blocks", file.path(WORK_DIR, "qtl_block_significance.csv"))
EQTL_FILE     <- get_arg("--eqtl",   file.path(WORK_DIR, "eqtl.csv"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("  QTL051-054 x trans-chromosomal FL/FS QTL\n")
cat("  Genetic effects (coding combinations + epistasis)\n")
cat("============================================================\n\n")

# ---- Step 1: eQTL analysis - trans-eQTL targets among FL/FS QTLs ----
cat("1. eQTL analysis: QTL051-054 -> trans-chromosomal FL/FS regulation...\n")

blocks <- fread(BLOCK_FILE)
eqtl <- fread(EQTL_FILE)

# Trans-eQTL from A07 genes
eqtl_a07_trans <- eqtl[grepl("Ghir_A07", `Phenotype ID`) & Type == "trans"]

# Parse SNP location
eqtl_a07_trans$snp_chr <- sub("_.*", "", eqtl_a07_trans$`QTL site`)
eqtl_a07_trans$snp_pos <- as.numeric(sub(".*_", "", eqtl_a07_trans$`QTL site`))

# Map SNP to blocks
map_to_block <- function(chr, pos, blk_df) {
  sapply(seq_along(chr), function(i) {
    bc <- chr[i]; bp <- pos[i]
    hits <- blk_df[grepl(bc, chr, fixed=TRUE) & block_start <= bp & block_end >= bp]
    if (nrow(hits) > 0) paste(unique(hits$QTL_id), collapse=";") else NA
  })
}
eqtl_a07_trans$target_QTL <- map_to_block(eqtl_a07_trans$snp_chr, eqtl_a07_trans$snp_pos, blocks)

mapped <- eqtl_a07_trans[!is.na(target_QTL)]

# Expand multi-QTL mappings
qtl_list <- strsplit(mapped$target_QTL, ";")
expanded <- data.table(
  gene = rep(mapped$`Phenotype ID`, lengths(qtl_list)),
  snp_chr = rep(mapped$snp_chr, lengths(qtl_list)),
  snp_pos = rep(mapped$snp_pos, lengths(qtl_list)),
  target_QTL = unlist(qtl_list),
  Stage = rep(mapped$Stage, lengths(qtl_list)),
  P_value = rep(mapped$`P-value`, lengths(qtl_list)),
  P_adjust = rep(mapped$`P-adjust`, lengths(qtl_list)),
  R2 = rep(mapped$`Adjusted R-square`, lengths(qtl_list))
)
expanded <- unique(expanded)

# FL/FS QTLs on non-A07 chromosomes
fl_fs_qtls <- unique(blocks[grepl("FL|FS", Sources) & !grepl("A07", chr), QTL_id])
fl_fs_expanded <- expanded[target_QTL %in% fl_fs_qtls]

# Summary by target QTL
eqtl_summary <- fl_fs_expanded[, .(
  n_genes = length(unique(gene)),
  n_snp_sites = length(unique(snp_pos)),
  stages = paste(sort(unique(Stage)), collapse=","),
  min_P = min(P_value),
  max_R2 = max(R2)
), by = .(target_QTL, snp_chr)]
eqtl_summary <- eqtl_summary[order(min_P)]

cat(sprintf("   QTL051-054 trans-eQTL targeting non-A07 FL/FS QTLs: %d records\n", nrow(fl_fs_expanded)))
cat(sprintf("   involving %d remote QTLs\n\n", length(unique(fl_fs_expanded$target_QTL))))
print(eqtl_summary)

fwrite(eqtl_summary, file.path(OUT_DIR, "trans_eQTL_FL_FS_summary.csv"))

# ---- Step 2: top QTL pairs, consensus coding combination analysis ----
cat("\n2. Consensus coding combination analysis...\n")

# Classify the coding state of an accession across all blocks of a QTL
classify_accession <- function(mat) {
  apply(mat, 1, function(x) {
    vals <- x[!is.na(x)]
    if (length(vals) == 0) return("NA")
    has2 <- any(vals == 2); has0 <- any(vals == 0)
    if (has2 && !has0) return("two")
    if (!has2 && has0) return("zero")
    if (has2 && has0) return("drag")
    return("NA")
  })
}

# QTL051-054 FL/FS traits
local_qtls <- list(
  QTL051 = list(traits = c("FL","FS")),
  QTL053 = list(traits = c("FS")),
  QTL054 = list(traits = c("FL","FS"))
)

# All QTLs with FL/FS consensus coding (discovered automatically)
all_consensus_folders <- list.dirs(CONSENSUS_DIR, recursive=FALSE, full.names=FALSE)
remote_qtls <- c()
for (f in all_consensus_folders) {
  files <- list.files(file.path(CONSENSUS_DIR, f), pattern="consensus_coding")
  if (any(grepl("^(FL|FS)_consensus", files))) {
    remote_qtls <- c(remote_qtls, f)
  }
}
# Exclude QTL051-054 themselves and QTL055 (also on A07)
remote_qtls <- setdiff(remote_qtls, c("QTL_051","QTL_052","QTL_053","QTL_054","QTL_055"))
cat(sprintf("   total remote FL/FS QTLs: %d\n", length(remote_qtls)))

read_consensus <- function(qtl, trait) {
  fn <- file.path(CONSENSUS_DIR, qtl, paste0(trait, "_consensus_coding.csv"))
  if (file.exists(fn)) {
    x <- fread(fn)
    blk_cols <- grep("^Block_", names(x), value=TRUE)
    x$class <- classify_accession(as.matrix(x[, ..blk_cols]))
    x$source_QTL <- qtl
    x$source_trait <- trait
    # Keep Accession, bio, phenotype, class
    pheno_name <- trait
    pheno_val <- if(pheno_name %in% names(x)) x[[pheno_name]] else NA_real_
    out <- data.table(Accession = x$Accession, bio = x$bio, class = x$class,
                      pheno = pheno_val)
    out
  } else NULL
}

# Read all local and remote consensus matrices
all_data <- list()
for (qtl_name in names(local_qtls)) {
  for (trait in local_qtls[[qtl_name]]$traits) {
    qtl_folder <- sub("0", "_0", sub("(\\d)(\\d)", "\\1_\\2", qtl_name))
    # QTL051 -> QTL_051
    qtl_folder <- paste0("QTL_", sub("QTL", "", qtl_name))
    dat <- read_consensus(qtl_folder, trait)
    if (!is.null(dat)) {
      key <- paste0(qtl_name, "_", trait)
      setnames(dat, c("class","pheno"), c(paste0("class_", key), paste0("pheno_", key)))
      all_data[[key]] <- dat
      cat(sprintf("   %s %s: n=%d\n", qtl_name, trait, nrow(dat)))
    }
  }
}

for (qtl_folder in remote_qtls) {
  # Generate QTL name: QTL_017 -> QTL017
  qtl_name <- gsub("_", "", qtl_folder)
  for (trait in c("FL", "FS")) {
    dat <- read_consensus(qtl_folder, trait)
    if (!is.null(dat)) {
      key <- paste0(qtl_name, "_", trait)
      setnames(dat, c("class","pheno"), c(paste0("class_", key), paste0("pheno_", key)))
      all_data[[key]] <- dat
      cat(sprintf("   %s %s: n=%d (remote)\n", qtl_name, trait, nrow(dat)))
    }
  }
}

# ---- Step 3: merge and combination-type analysis ----
cat("\n3. Merging and combination analysis...\n")

# Merge all data into a single table
merged <- NULL
for (nm in names(all_data)) {
  dat <- all_data[[nm]]
  if (is.null(merged)) {
    merged <- dat
  } else {
    merged <- merge(merged, dat, by=c("Accession","bio"), all=TRUE)
  }
}
cat(sprintf("   total accessions: %d\n", nrow(merged)))

# Analyse all local x remote combinations
all_results <- list()
all_pheno <- list()

local_keys <- grep("^class_QTL05[1-4]_", names(merged), value=TRUE)
remote_keys <- grep("^class_QTL", names(merged), value=TRUE)
remote_keys <- setdiff(remote_keys, local_keys)

target_pops <- c("Group-A", "Group-B", "Group-C")

for (lk in local_keys) {
  for (rk in remote_keys) {
    # coding classes
    sub <- merged[!is.na(get(lk)) & !is.na(get(rk))]
    sub <- sub[get(lk) %in% c("two","zero") & get(rk) %in% c("two","zero")]

    if (nrow(sub) < 30) next

    sub$combo <- paste(sub[[lk]], sub[[rk]], sep="_")

    # combination naming
    sub$combo_type <- with(sub, case_when(
      get(lk) == "two" & get(rk) == "two" ~ "Both_Superior",
      get(lk) == "zero" & get(rk) == "zero" ~ "Both_Inferior",
      get(lk) == "two" & get(rk) == "zero" ~ "Local_Sup_Remote_Inf",
      get(lk) == "zero" & get(rk) == "two" ~ "Local_Inf_Remote_Sup",
      TRUE ~ "Other"
    ))

    # local trait name
    local_trait <- sub("class_", "", lk)  # e.g., "QTL051_FL"
    local_short <- sub("QTL\\d+_", "", lk)  # e.g., "FL" or "FS"

    # phenotypic validation
    pheno_col <- paste0("pheno_", local_trait)
    if (pheno_col %in% names(sub)) {
      sub_pheno <- sub[!is.na(get(pheno_col))]

      if (nrow(sub_pheno) >= 20) {
        smry <- sub_pheno[, .(
          n = .N,
          pheno_mean = mean(get(pheno_col), na.rm=TRUE),
          pheno_sd = sd(get(pheno_col), na.rm=TRUE)
        ), by = combo_type]

        dS <- sub_pheno[combo_type == "Both_Superior"]
        dI <- sub_pheno[combo_type == "Both_Inferior"]

        pval <- NA_real_; epistasis <- NA_real_
        if (nrow(dS) >= 3 && nrow(dI) >= 3) {
          tt <- t.test(dS[[pheno_col]], dI[[pheno_col]])
          pval <- tt$p.value
        }

        # epistasis test: two-way ANOVA
        sub_pheno$G_local <- ifelse(sub_pheno[[lk]] == "two", 1, 0)
        sub_pheno$G_remote <- ifelse(sub_pheno[[rk]] == "two", 1, 0)
        sub_pheno$y <- sub_pheno[[pheno_col]]
        aov_fit <- tryCatch(
          summary(aov(y ~ G_local * G_remote, data=sub_pheno)),
          error = function(e) NULL
        )
        aov_interact_p <- NA_real_
        if (!is.null(aov_fit)) {
          aov_interact_p <- aov_fit[[1]]$`Pr(>F)`[3]
        }

        # record
        pair_name <- paste(lk, rk, sep=" x ")
        all_pheno[[pair_name]] <- data.table(
          local_QTL = lk, remote_QTL = rk,
          n_total = nrow(sub_pheno),
          n_Both_Sup = nrow(dS),
          n_Both_Inf = nrow(dI),
          mean_Sup = if(nrow(dS)>0) mean(dS[[pheno_col]], na.rm=TRUE) else NA,
          mean_Inf = if(nrow(dI)>0) mean(dI[[pheno_col]], na.rm=TRUE) else NA,
          delta = if(nrow(dS)>0 && nrow(dI)>0)
            mean(dS[[pheno_col]],na.rm=TRUE) - mean(dI[[pheno_col]],na.rm=TRUE) else NA,
          ttest_P = pval,
          anova_interact_P = aov_interact_p
        )

        # print key results
        cat(sprintf("\n   %s x %s [%s]:\n", lk, rk, local_trait))
        cat(sprintf("      Both Superior: n=%d, %.3f +/- %.3f\n",
          nrow(dS), if(nrow(dS)>0) mean(dS[[pheno_col]],na.rm=TRUE) else NA,
          if(nrow(dS)>0) sd(dS[[pheno_col]],na.rm=TRUE) else NA))
        cat(sprintf("      Both Inferior: n=%d, %.3f +/- %.3f\n",
          nrow(dI), if(nrow(dI)>0) mean(dI[[pheno_col]],na.rm=TRUE) else NA,
          if(nrow(dI)>0) sd(dI[[pheno_col]],na.rm=TRUE) else NA))
        cat(sprintf("      t-test P = %.4f | ANOVA interaction P = %.4f\n", pval, aov_interact_p))

        # ---- per-population analysis ----
        for (pop in target_pops) {
          pop_sub <- sub_pheno[bio == pop]
          if (nrow(pop_sub) < 15) next
          pop_dS <- pop_sub[combo_type == "Both_Superior"]
          pop_dI <- pop_sub[combo_type == "Both_Inferior"]

          pop_pval <- NA_real_
          if (nrow(pop_dS) >= 3 && nrow(pop_dI) >= 3) {
            pop_tt <- t.test(pop_dS[[pheno_col]], pop_dI[[pheno_col]])
            pop_pval <- pop_tt$p.value
          }

          cat(sprintf("      [%s] S-S: n=%d, %.2f +/- %.2f | I-I: n=%d, %.2f +/- %.2f | P=%.4f\n",
            pop,
            nrow(pop_dS), if(nrow(pop_dS)>0) mean(pop_dS[[pheno_col]],na.rm=TRUE) else NA,
            if(nrow(pop_dS)>0) sd(pop_dS[[pheno_col]],na.rm=TRUE) else NA,
            nrow(pop_dI), if(nrow(pop_dI)>0) mean(pop_dI[[pheno_col]],na.rm=TRUE) else NA,
            if(nrow(pop_dI)>0) sd(pop_dI[[pheno_col]],na.rm=TRUE) else NA,
            pop_pval))
        }

        # store the detailed classification table
        smry$local_QTL <- lk; smry$remote_QTL <- rk; smry$trait <- local_trait
        all_results[[pair_name]] <- smry
      }
    }
  }
}

# ---- Step 4: output ----
cat("\n\n4. Writing output...\n")

all_pheno_dt <- rbindlist(all_pheno, fill=TRUE)
all_pheno_dt <- all_pheno_dt[order(anova_interact_P)]
fwrite(all_pheno_dt, file.path(OUT_DIR, "cross_QTL_epistasis_summary.csv"))
cat(sprintf("   epistasis summary: cross_QTL_epistasis_summary.csv (%d rows)\n", nrow(all_pheno_dt)))

all_res_dt <- rbindlist(all_results, fill=TRUE)
fwrite(all_res_dt, file.path(OUT_DIR, "cross_QTL_combo_pheno.csv"))
cat(sprintf("   combination phenotypes: cross_QTL_combo_pheno.csv (%d rows)\n", nrow(all_res_dt)))

# full merged table
fwrite(merged, file.path(OUT_DIR, "cross_QTL_all_merged.csv"))
cat(sprintf("   full merged table: cross_QTL_all_merged.csv\n"))

# ---- Step 5: summary of key findings ----
cat("\n============================================================\n")
cat("  Key findings\n")
cat("============================================================\n")

if (nrow(all_pheno_dt) > 0) {
  cat("\n  Significant cross-chromosomal QTL interactions (ANOVA interaction P < 0.05):\n\n")
  sig_pairs <- all_pheno_dt[anova_interact_P < 0.05]
  if (nrow(sig_pairs) > 0) {
    for (i in seq_len(min(10, nrow(sig_pairs)))) {
      sp <- sig_pairs[i]
      cat(sprintf("   %s x %s\n", sp$local_QTL, sp$remote_QTL))
      cat(sprintf("      Both_Sup=%.2f  Both_Inf=%.2f  delta=%.2f\n",
        sp$mean_Sup, sp$mean_Inf, sp$delta))
      cat(sprintf("      ANOVA interaction P = %.2e  t-test P = %.4f\n",
        sp$anova_interact_P, sp$ttest_P))
    }
  } else {
    cat("   (no significant epistatic signal)\n")
  }

  cat("\n  Top 10 QTL pairs (by effect size):\n\n")
  top10 <- all_pheno_dt[order(-abs(delta))][1:min(10, nrow(all_pheno_dt))]
  print(top10[, .(local_QTL, remote_QTL, n_total, mean_Sup, mean_Inf, delta,
    ttest_P, anova_interact_P)])
}

cat("\nDone.\n")
