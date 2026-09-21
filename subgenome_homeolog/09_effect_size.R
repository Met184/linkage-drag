#!/usr/bin/env Rscript
# 09_effect_size.R
# Effect size of the original block and of its homeologous interval, and their ratio.
# Usage: Rscript 09_effect_size.R [--home DIR] [--detail DIR] [--out DIR]

suppressMessages({
  library(dplyr)
})

# ---- CONFIG -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1L) args[i + 1L] else default
}

HOME    <- get_arg("--home", ".")
DETAIL  <- get_arg("--detail", file.path(HOME, "detail"))
OUT_DIR <- get_arg("--out", file.path(HOME, "effect_size"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Blocks to analyse. Each case names the detail table and the trait it holds.
cases <- list(
  list(block = "Block_123", qtl = "QTL_051", trait = "FS"),
  list(block = "Block_146", qtl = "QTL_054", trait = "FS"),
  list(block = "Block_155", qtl = "QTL_054", trait = "FL"),
  list(block = "Block_332", qtl = "QTL_127", trait = "FD")
)

wtest <- function(a, b) {
  if (length(a) < 3 || length(b) < 3) return(NA_real_)
  tryCatch(wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA_real_)
}

res <- list()
for (case in cases) {
  d <- read.csv(file.path(DETAIL, paste0("detail_", case$block, ".csv")),
                stringsAsFactors = FALSE, check.names = FALSE)
  v <- case$trait
  bios <- c(sort(unique(d$bio)), "ALL")
  for (b in bios) {
    dd <- if (b == "ALL") d else d[d$bio == b, ]

    # ---- original effect: superior (2) vs inferior (0) ----
    xS <- dd[[v]][dd$orig_code == 2]; xI <- dd[[v]][dd$orig_code == 0]
    mS <- mean(xS, na.rm = TRUE); mI <- mean(xI, na.rm = TRUE)
    nS <- length(xS); nI <- length(xI)
    d_orig  <- mS - mI
    p_orig  <- wtest(xS, xI)

    # ---- homeolog effect (conditional on an inferior original block) ----
    d0 <- dd[dd$orig_code == 0, ]
    yS <- d0[[v]][d0$homeo_code == 2]; yI <- d0[[v]][d0$homeo_code == 0]
    mhS <- mean(yS, na.rm = TRUE); mhI <- mean(yI, na.rm = TRUE)
    nhS <- length(yS); nhI <- length(yI)
    d_homeo  <- mhS - mhI
    p_homeo  <- wtest(yS, yI)

    # ---- ratio ----
    ratio_signed <- if (is.na(d_homeo) || d_homeo == 0) NA_real_ else d_orig / d_homeo
    ratio_abs    <- if (is.na(d_homeo) || d_homeo == 0) NA_real_ else abs(d_orig) / abs(d_homeo)

    res[[length(res) + 1]] <- data.frame(
      block = case$block, qtl = case$qtl, trait = case$trait, bio = b,
      orig_n_S = nS, orig_mean_S = mS, orig_n_I = nI, orig_mean_I = mI,
      orig_effect = d_orig, orig_wilcox_p = p_orig,
      homeo_n_S = nhS, homeo_mean_S = mhS, homeo_n_I = nhI, homeo_mean_I = mhI,
      homeo_effect = d_homeo, homeo_wilcox_p = p_homeo,
      ratio_signed = ratio_signed, ratio_abs = ratio_abs,
      stringsAsFactors = FALSE)
  }
}

out <- do.call(rbind, res)
out[, -c(1:4)] <- lapply(out[, -c(1:4)], function(x) signif(x, 4))
write.csv(out, file.path(OUT_DIR, "effect_size_results.csv"), row.names = FALSE)

cat("=== effect size results ===\n")
print(out, row.names = FALSE)
cat("\nwritten to:", file.path(OUT_DIR, "effect_size_results.csv"), "\n")
