#!/usr/bin/env Rscript
# 02_weighted_R_load.R
# Per-accession unweighted and weighted antagonistic (R) load.
# Usage: Rscript 02_weighted_R_load.R [--work DIR] [--out DIR]

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(rstatix))
suppressMessages(library(multcompView))

# ---- CONFIG -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1L) args[i + 1L] else default
}

WORK       <- get_arg("--work", ".")
ENCODE     <- get_arg("--encode",     file.path(WORK, "encode"))
DIR_R      <- get_arg("--antagonism", file.path(WORK, "antagonism"))
DIR_BIO    <- get_arg("--pheno",      file.path(WORK, "bio"))
META       <- get_arg("--meta",       file.path(WORK, "meta"))
OUT        <- get_arg("--out",        file.path(WORK, "weighted"))
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# Antagonism classes used to split the R load by QTL antagonism rate
CAT_HIGH <- "High"
CAT_LOW  <- "Low"
CAT_NONE <- "None"

# =====================================================================
# Step 1: Load phenotype data
# =====================================================================
bio_files <- list.files(DIR_BIO, pattern="\\.csv$", full.names=TRUE)
bio_all <- rbindlist(lapply(bio_files, function(f) {
  trait <- sub("\\.csv$", "", basename(f))
  d <- fread(f)
  setnames(d, c("Accession","Population","Phenotype"))
  d[, Trait := trait]
  d[!is.na(Phenotype)]
}))

# =====================================================================
# Step 2: Load SIR + R-type
# =====================================================================
r_qtl <- fread(file.path(DIR_R, "antagonism_by_QTL.csv"))
r_qtl[, R_Group := fcase(Mean_R >= 20, CAT_HIGH, Mean_R > 0 & Mean_R < 20, CAT_LOW, default = CAT_NONE)]
qtl_r_group <- setNames(r_qtl$R_Group, r_qtl$QTL)

files <- list.files(ENCODE, pattern="_SIR\\.csv$", full=FALSE)
sir_hash <- list()
for (f in files) {
  base <- sub("_SIR\\.csv$", "", f)
  parts <- strsplit(base, "_")[[1]]
  qnum <- parts[2]; qtl <- paste0("QTL", qnum)
  trait <- if(length(parts) > 3) paste(parts[3:length(parts)], collapse="_") else parts[3]
  d <- fread(file.path(ENCODE, f))
  accs <- as.character(d[[1]]); codes <- as.character(d[[2]])
  sir_hash[[qtl]][[trait]] <- setNames(codes, accs)
}
all_accs <- unique(unlist(lapply(sir_hash, function(q) names(q[[1]]))))

# =====================================================================
# Step 3: Compute |d| for ALL QTL-trait pairs
# =====================================================================
cat("Computing effect sizes for all QTL-traits...\n")
qtl_traits_all <- rbindlist(lapply(names(sir_hash), function(qtl) {
  data.table(QTL=qtl, Trait=names(sir_hash[[qtl]]))
}))
cat(sprintf("  Total: %d QTL-trait pairs\n", nrow(qtl_traits_all)))

compute_d <- function(qtl, trait) {
  codes <- sir_hash[[qtl]][[trait]]
  if (is.null(codes)) return(NA_real_)
  bio_trait <- bio_all[Trait == trait]
  if (nrow(bio_trait) == 0) return(NA_real_)
  sir_dt <- data.table(Accession=names(codes), SIR=as.character(codes))
  merged <- merge(sir_dt, bio_trait, by="Accession")
  s_accs <- merged[SIR == "S"]; i_accs <- merged[SIR == "I"]
  if (nrow(s_accs) < 5 || nrow(i_accs) < 5) return(NA_real_)
  mean_S <- mean(s_accs$Phenotype); mean_I <- mean(i_accs$Phenotype)
  sd_S <- sd(s_accs$Phenotype); sd_I <- sd(i_accs$Phenotype)
  sd_pooled <- sqrt(((nrow(s_accs)-1)*sd_S^2 + (nrow(i_accs)-1)*sd_I^2) / (nrow(s_accs)+nrow(i_accs)-2))
  abs((mean_S - mean_I) / sd_pooled)
}

qtl_traits_all[, d_abs := mapply(compute_d, QTL, Trait)]
cat(sprintf("  Valid d: %d / NA: %d\n\n", sum(!is.na(qtl_traits_all$d_abs)), sum(is.na(qtl_traits_all$d_abs))))

# =====================================================================
# Step 4: Per-accession R counting - unweighted + weighted
# =====================================================================
cat("Counting R (unweighted + weighted)...\n")

# For weighted: normalize d_abs so mean weight = 1
valid_d <- qtl_traits_all[!is.na(d_abs)]
valid_d[, d_norm := d_abs / mean(d_abs)]
# Build lookup: QTL|Trait -> normalized d
d_lookup <- setNames(valid_d$d_norm, paste0(valid_d$QTL, "|", valid_d$Trait))

init <- function() list(R_high_unwt=0, R_low_unwt=0, R_high_wt=0, R_low_wt=0,
                        n_valid_high=0L, n_valid_low=0L, n_valid_all=0L)
acc <- setNames(lapply(all_accs, function(a) init()), all_accs)

for (qtl in names(sir_hash)) {
  grp <- qtl_r_group[qtl]; if (is.na(grp)) grp <- CAT_NONE
  for (trait in names(sir_hash[[qtl]])) {
    codes <- sir_hash[[qtl]][[trait]]
    r_accs <- names(codes)[codes == "R"]
    w <- d_lookup[paste0(qtl, "|", trait)]
    wt <- if (is.na(w)) 0 else w
    if (grp == CAT_HIGH) {
      valid_accs <- names(codes)[codes %in% c("S","I","R")]
      for (a in valid_accs) {
        acc[[a]]$n_valid_high <- acc[[a]]$n_valid_high + 1L
        acc[[a]]$n_valid_all  <- acc[[a]]$n_valid_all + 1L
      }
      for (a in r_accs) {
        acc[[a]]$R_high_unwt <- acc[[a]]$R_high_unwt + 1L
        acc[[a]]$R_high_wt   <- acc[[a]]$R_high_wt + wt
      }
    } else if (grp == CAT_LOW) {
      valid_accs <- names(codes)[codes %in% c("S","I","R")]
      for (a in valid_accs) {
        acc[[a]]$n_valid_low <- acc[[a]]$n_valid_low + 1L
        acc[[a]]$n_valid_all  <- acc[[a]]$n_valid_all + 1L
      }
      for (a in r_accs) {
        acc[[a]]$R_low_unwt <- acc[[a]]$R_low_unwt + 1L
        acc[[a]]$R_low_wt   <- acc[[a]]$R_low_wt + wt
      }
    } else {
      valid_accs <- names(codes)[codes %in% c("S","I","R")]
      for (a in valid_accs) {
        acc[[a]]$n_valid_all  <- acc[[a]]$n_valid_all + 1L
      }
    }
  }
}

acc_df <- rbindlist(lapply(names(acc), function(a) {
  v <- acc[[a]]
  data.table(Accession=a,
    R_high_unwt=v$R_high_unwt, R_low_unwt=v$R_low_unwt,
    R_total_unwt=v$R_high_unwt+v$R_low_unwt,
    R_high_wt=v$R_high_wt, R_low_wt=v$R_low_wt,
    R_total_wt=v$R_high_wt+v$R_low_wt,
    n_valid_high=v$n_valid_high, n_valid_low=v$n_valid_low, n_valid_all=v$n_valid_all)
}))
# Rates
acc_df[, `:=`(
  R_total_unwt_pct = ifelse(n_valid_all>0, R_total_unwt/n_valid_all*100, NA),
  R_total_wt_pct   = ifelse(n_valid_all>0, R_total_wt/n_valid_all*100, NA),
  R_high_unwt_pct  = ifelse(n_valid_high>0, R_high_unwt/n_valid_high*100, NA),
  R_high_wt_pct    = ifelse(n_valid_high>0, R_high_wt/n_valid_high*100, NA),
  R_low_unwt_pct   = ifelse(n_valid_low>0, R_low_unwt/n_valid_low*100, NA),
  R_low_wt_pct     = ifelse(n_valid_low>0, R_low_wt/n_valid_low*100, NA)
)]

cat(sprintf("Mean R_total: unweighted=%.2f  weighted=%.2f\n",
  mean(acc_df$R_total_unwt), mean(acc_df$R_total_wt)))

# =====================================================================
# Step 5: Meta + era
# =====================================================================
meta <- fread(file.path(META, "yangbenALL.csv"))
meta[, Accession := as.character(Accession)]
meta <- meta[, .(Accession, Type, Date)]

acc_all <- merge(acc_df, meta, by="Accession", all.x=TRUE)
acc_all[, Group := fcase(
  Type == "Semi-wild", "Semi-wild",
  Date < 1950, "<1950",
  Date >= 1950 & Date < 1960, "1950-1959",
  Date >= 1960 & Date < 1970, "1960-1969",
  Date >= 1970 & Date < 1980, "1970-1979",
  Date >= 1980 & Date < 1990, "1980-1989",
  Date >= 1990 & Date < 2000, "1990-1999",
  Date >= 2000 & Date < 2010, "2000-2009",
  Date >= 2010, "2010-2017",
  default = NA_character_
)]
groups <- c("Semi-wild","<1950","1950-1959","1960-1969","1970-1979",
            "1980-1989","1990-1999","2000-2009","2010-2017")
acc_plot <- acc_all[Group %in% groups]
acc_plot[, Group := factor(Group, levels=groups)]

# =====================================================================
# Step 6: Comparison
# =====================================================================
fmt_p <- function(p) {
  if (p < 0.0001) return("p < 0.0001")
  if (p < 0.001)  return(sprintf("p = %.4f", p))
  if (p < 0.01)   return(sprintf("p = %.4f", p))
  if (p < 0.05)   return(sprintf("p = %.4f", p))
  return(sprintf("p = %.3f (ns)", p))
}

cat("\n========== R COMPARISON: Unweighted vs Weighted ==========\n")
cat(sprintf("%-14s %12s %12s %8s %12s %12s %8s %12s %12s %8s\n",
  "Era", "R_high_u", "R_high_w", "delta", "R_low_u", "R_low_w", "delta", "R_total_u", "R_total_w", "delta"))
cat(strrep("-", 100), "\n")
for(g in groups) {
  s <- acc_plot[Group == g]
  rhu <- mean(s$R_high_unwt); rhw <- mean(s$R_high_wt)
  rlu <- mean(s$R_low_unwt); rlw <- mean(s$R_low_wt)
  rtu <- mean(s$R_total_unwt); rtw <- mean(s$R_total_wt)
  cat(sprintf("%-14s %11.2f %11.2f %+7.2f %11.2f %11.2f %+7.2f %11.2f %11.2f %+7.2f\n",
    g, rhu, rhw, rhw-rhu, rlu, rlw, rlw-rlu, rtu, rtw, rtw-rtu))
}

cat("\n========== KW ==========\n")
for (yv in c("R_total_unwt","R_total_wt","R_high_unwt","R_high_wt","R_low_unwt","R_low_wt")) {
  kw <- kruskal.test(as.formula(paste0(yv, " ~ Group")), data=acc_plot)
  cat(sprintf("%s: %s\n", yv, fmt_p(kw$p.value)))
}

# =====================================================================
# Step 7: Plot
# =====================================================================
get_letters <- function(df, y_var) {
  f <- as.formula(paste0(y_var, " ~ Group"))
  dunn <- dunn_test(df, f, p.adjust.method = "bonferroni")
  setDT(dunn)
  grps <- as.character(unique(df$Group))
  pmat <- matrix(1, nrow=length(grps), ncol=length(grps), dimnames=list(grps, grps))
  for (i in 1:nrow(dunn)) {
    pmat[dunn$group1[i], dunn$group2[i]] <- dunn$p.adj[i]
    pmat[dunn$group2[i], dunn$group1[i]] <- dunn$p.adj[i]
  }
  lets <- multcompLetters(pmat, threshold=0.05)$Letters
  data.table(Group=names(lets), Letter=lets)
}

white_theme <- theme_bw() + theme(
  panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
  panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),
  axis.title=element_text(size=12, face="bold"), axis.text=element_text(size=10),
  axis.text.x=element_text(size=8, angle=40, hjust=1),
  plot.title=element_text(size=12, face="bold", hjust=0.5), legend.position="none")

group_cols <- c("Semi-wild"="#8DD3C7","<1950"="#FFFFB3","1950-1959"="#BEBADA",
  "1960-1969"="#FB8072","1970-1979"="#80B1D3","1980-1989"="#FDB462",
  "1990-1999"="#B3DE69","2000-2009"="#FCCDE5","2010-2017"="#D9D9D9")

n_labs <- acc_plot[, .N, by=Group][order(Group)]
n_vec  <- setNames(paste0(n_labs$Group, "\n(n=", n_labs$N, ")"), n_labs$Group)

make_R_cmp_plot <- function(y_unwt, y_wt, title_text, y_label="R count") {
  d1 <- acc_plot[, .(Group, Value=get(y_unwt))]
  d1[, Method := "Unweighted"]
  d2 <- acc_plot[, .(Group, Value=get(y_wt))]
  d2[, Method := "Weighted (|d|)"]
  dm <- rbind(d1, d2)
  dm[, Method := factor(Method, levels=c("Unweighted","Weighted (|d|)"))]

  m_list <- list()
  for (mth in unique(dm$Method)) {
    sub <- dm[Method == mth]
    sub[, Group := factor(Group, levels=groups)]
    mm <- sub[, .(Mean=mean(Value, na.rm=TRUE), yMax=max(Value, na.rm=TRUE)), by=Group]
    ll <- get_letters(sub, "Value")
    mm <- merge(mm, ll, by="Group")
    yrange <- diff(range(c(0, mm$yMax)))
    if (yrange == 0 || is.na(yrange)) yrange <- 1
    mm[, LabelY := yMax + yrange * 0.06]
    mm[, Method := mth]
    m_list[[mth]] <- mm
  }
  m_all <- rbindlist(m_list)

  ggplot(dm, aes(x=Group, y=Value, fill=Group)) +
    geom_boxplot(outlier.size=0.3, alpha=0.85, width=0.6) +
    geom_label(data=m_all, aes(x=Group, y=LabelY, label=sprintf("%.1f %s", Mean, Letter)),
               fill="white", alpha=0.85, size=2.3, fontface="bold",
               vjust=-0.5, label.padding=unit(0.1, "lines")) +
    facet_wrap(~Method, ncol=1, scales="free_y") +
    scale_fill_manual(values=group_cols) + scale_x_discrete(labels=n_vec) +
    labs(x="", y=y_label, title=title_text) + white_theme
}

p_rh <- make_R_cmp_plot("R_high_unwt", "R_high_wt", "R_high (high-antagonism QTLs)")
p_rl <- make_R_cmp_plot("R_low_unwt", "R_low_wt", "R_low (low-antagonism QTLs)")
p_rt <- make_R_cmp_plot("R_total_unwt", "R_total_wt", "R_total")

ggsave(file.path(OUT, "R_high_weighted_vs_unweighted.pdf"), p_rh, width=10, height=7)
ggsave(file.path(OUT, "R_high_weighted_vs_unweighted.png"), p_rh, width=10, height=7, dpi=300)
ggsave(file.path(OUT, "R_low_weighted_vs_unweighted.pdf"), p_rl, width=10, height=7)
ggsave(file.path(OUT, "R_low_weighted_vs_unweighted.png"), p_rl, width=10, height=7, dpi=300)
ggsave(file.path(OUT, "R_total_weighted_vs_unweighted.pdf"), p_rt, width=10, height=7)
ggsave(file.path(OUT, "R_total_weighted_vs_unweighted.png"), p_rt, width=10, height=7, dpi=300)

# --- Rate plots ---
p_rh_rate <- make_R_cmp_plot("R_high_unwt_pct", "R_high_wt_pct", "R_high rate (high-antagonism QTLs)", y_label="R rate (%)")
p_rl_rate <- make_R_cmp_plot("R_low_unwt_pct", "R_low_wt_pct", "R_low rate (low-antagonism QTLs)", y_label="R rate (%)")
p_rt_rate <- make_R_cmp_plot("R_total_unwt_pct", "R_total_wt_pct", "R_total rate", y_label="R rate (%)")

ggsave(file.path(OUT, "R_high_rate_weighted_vs_unweighted.pdf"), p_rh_rate, width=10, height=7)
ggsave(file.path(OUT, "R_high_rate_weighted_vs_unweighted.png"), p_rh_rate, width=10, height=7, dpi=300)
ggsave(file.path(OUT, "R_low_rate_weighted_vs_unweighted.pdf"), p_rl_rate, width=10, height=7)
ggsave(file.path(OUT, "R_low_rate_weighted_vs_unweighted.png"), p_rl_rate, width=10, height=7, dpi=300)
ggsave(file.path(OUT, "R_total_rate_weighted_vs_unweighted.pdf"), p_rt_rate, width=10, height=7)
ggsave(file.path(OUT, "R_total_rate_weighted_vs_unweighted.png"), p_rt_rate, width=10, height=7, dpi=300)

# --- KW for rates ---
cat("\n========== KW (Rates) ==========\n")
for (yv in c("R_total_unwt_pct","R_total_wt_pct","R_high_unwt_pct","R_high_wt_pct","R_low_unwt_pct","R_low_wt_pct")) {
  kw <- kruskal.test(as.formula(paste0(yv, " ~ Group")), data=acc_plot)
  cat(sprintf("%s: %s\n", yv, fmt_p(kw$p.value)))
}

# --- Era comparison for rates ---
cat("\n========== R RATE COMPARISON: Unweighted vs Weighted ==========\n")
cat(sprintf("%-14s %12s %12s %8s %12s %12s %8s %12s %12s %8s\n",
  "Era", "R_high_u%", "R_high_w%", "delta", "R_low_u%", "R_low_w%", "delta", "R_total_u%", "R_total_w%", "delta"))
cat(strrep("-", 100), "\n")
for(g in groups) {
  s <- acc_plot[Group == g]
  rhu <- mean(s$R_high_unwt_pct, na.rm=TRUE); rhw <- mean(s$R_high_wt_pct, na.rm=TRUE)
  rlu <- mean(s$R_low_unwt_pct, na.rm=TRUE); rlw <- mean(s$R_low_wt_pct, na.rm=TRUE)
  rtu <- mean(s$R_total_unwt_pct, na.rm=TRUE); rtw <- mean(s$R_total_wt_pct, na.rm=TRUE)
  cat(sprintf("%-14s %11.2f %11.2f %+7.2f %11.2f %11.2f %+7.2f %11.2f %11.2f %+7.2f\n",
    g, rhu, rhw, rhw-rhu, rlu, rlw, rlw-rlu, rtu, rtw, rtw-rtu))
}

# --- Save summaries (absolute + rate) ---
summ_abs <- acc_plot[, .(N=.N,
  R_total_u=round(mean(R_total_unwt),2), R_total_w=round(mean(R_total_wt),2),
  R_high_u=round(mean(R_high_unwt),2),  R_high_w=round(mean(R_high_wt),2),
  R_low_u=round(mean(R_low_unwt),2),    R_low_w=round(mean(R_low_wt),2)
), by=Group][order(Group)]

summ_rate <- acc_plot[, .(N=.N,
  R_total_u_pct=round(mean(R_total_unwt_pct, na.rm=TRUE),1),
  R_total_w_pct=round(mean(R_total_wt_pct, na.rm=TRUE),1),
  R_high_u_pct=round(mean(R_high_unwt_pct, na.rm=TRUE),1),
  R_high_w_pct=round(mean(R_high_wt_pct, na.rm=TRUE),1),
  R_low_u_pct=round(mean(R_low_unwt_pct, na.rm=TRUE),1),
  R_low_w_pct=round(mean(R_low_wt_pct, na.rm=TRUE),1)
), by=Group][order(Group)]

cat("\n\n========== SUMMARY (Absolute) ==========\n")
print(summ_abs)
cat("\n\n========== SUMMARY (Rate %) ==========\n")
print(summ_rate)

fwrite(summ_abs,  file.path(OUT, "R_weighted_vs_unweighted_by_period_counts.csv"), bom=TRUE)
fwrite(summ_rate, file.path(OUT, "R_weighted_vs_unweighted_by_period_rates.csv"), bom=TRUE)
fwrite(acc_all, file.path(OUT, "R_per_accession_weighted.csv"), bom=TRUE)
fwrite(qtl_traits_all, file.path(OUT, "all_QTL_trait_effect_sizes.csv"), bom=TRUE)

cat("\n=== Done ===\n")
