#!/usr/bin/env Rscript
# 01_weighted_drag_load.R
# Per-accession unweighted and weighted QTL-internal and QTL-between drag load.
# Usage: Rscript 01_weighted_drag_load.R [--work DIR] [--out DIR]

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
DIR_INT    <- get_arg("--internal",   file.path(WORK, "internal_drag"))
DIR_BET    <- get_arg("--between",    file.path(WORK, "between_drag"))
DIR_R      <- get_arg("--antagonism", file.path(WORK, "antagonism"))
DIR_BIO    <- get_arg("--pheno",      file.path(WORK, "bio"))
META       <- get_arg("--meta",       file.path(WORK, "meta"))
OUT        <- get_arg("--out",        file.path(WORK, "weighted"))
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# Antagonism classes used to split the R load by QTL antagonism rate
CAT_HIGH <- "High"
CAT_LOW  <- "Low"
CAT_NONE <- "None"
PURE_DRAG <- c("S-I", "I-S")

cat("============================================================\n")
cat("WEIGHTED LINKAGE DRAG ANALYSIS\n")
cat("Weights: product of haplotype effect sizes (Cohen's d)\n")
cat("============================================================\n\n")

# =====================================================================
# Step 1: Load bio (phenotype) data for all traits
# =====================================================================
cat("Loading phenotype data...\n")
bio_files <- list.files(DIR_BIO, pattern="\\.csv$", full.names=TRUE)
bio_all <- rbindlist(lapply(bio_files, function(f) {
  trait <- sub("\\.csv$", "", basename(f))
  d <- fread(f)
  setnames(d, c("Accession","Population","Phenotype"))
  d[, Trait := trait]
  d[!is.na(Phenotype)]
}))
cat(sprintf("  Loaded %d traits, %d records\n",
  uniqueN(bio_all$Trait), nrow(bio_all)))

# =====================================================================
# Step 2: Load SIR + filter pairs (>20% pure drag)
# =====================================================================
cat("Loading SIR encoding...\n")
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

# Internal pairs with pure drag rate
int_pairs_raw <- fread(file.path(DIR_INT, "QTL_trait_pairs.csv"), header=FALSE, skip=1)
int_pairs_all <- int_pairs_raw[, .(QTL=as.character(V1), Trait1=as.character(V3), Trait2=as.character(V4))]
int_drag_rates <- sapply(1:nrow(int_pairs_all), function(i) {
  qtl <- int_pairs_all$QTL[i]; t1 <- int_pairs_all$Trait1[i]; t2 <- int_pairs_all$Trait2[i]
  s1 <- sir_hash[[qtl]][[t1]]; s2 <- sir_hash[[qtl]][[t2]]
  if (is.null(s1) || is.null(s2)) return(NA_real_)
  common <- intersect(names(s1), names(s2))
  combos <- paste0(s1[common], "-", s2[common])
  valid <- nchar(s1[common]) == 1 & nchar(s2[common]) == 1
  if (sum(valid) == 0) return(0)
  sum(combos[valid] %in% PURE_DRAG) / sum(valid) * 100
})
int_pairs_all[, PureDragPct := round(int_drag_rates, 1)]
int_pairs <- int_pairs_all[PureDragPct > 20]

# Between pairs
bet_pairs_raw <- fread(file.path(DIR_BET, "QTL_between_drag_details.csv"))
bet_pairs_all <- bet_pairs_raw[, .(QTL1=as.character(bet_pairs_raw[[1]]), QTL2=as.character(bet_pairs_raw[[2]]),
                                Trait1=as.character(bet_pairs_raw[[3]]), Trait2=as.character(bet_pairs_raw[[4]]),
                                PureDragPct=as.numeric(bet_pairs_raw[[16]]))]
bet_pairs <- bet_pairs_all[PureDragPct > 20]

cat(sprintf("Filtered: internal=%d/%d  between=%d/%d\n\n",
  nrow(int_pairs), nrow(int_pairs_all), nrow(bet_pairs), nrow(bet_pairs_all)))

# =====================================================================
# Step 3: Compute effect sizes (Cohen's d) for each QTL-trait
# =====================================================================
cat("Computing haplotype effect sizes (Cohen's d)...\n")

# Build unique list of QTL-trait combos across all drag pairs
qtl_trait_pairs <- unique(rbind(
  rbindlist(list(
    int_pairs[, .(QTL=QTL, Trait=Trait1)],
    int_pairs[, .(QTL=QTL, Trait=Trait2)]
  )),
  rbindlist(list(
    bet_pairs[, .(QTL=QTL1, Trait=Trait1)],
    bet_pairs[, .(QTL=QTL1, Trait=Trait2)],
    bet_pairs[, .(QTL=QTL2, Trait=Trait1)],
    bet_pairs[, .(QTL=QTL2, Trait=Trait2)]
  ))
))
qtl_trait_pairs <- unique(qtl_trait_pairs)
cat(sprintf("  %d unique QTL-trait combinations to evaluate\n", nrow(qtl_trait_pairs)))

# Compute d for each QTL-trait
compute_d <- function(qtl, trait) {
  codes <- sir_hash[[qtl]][[trait]]
  if (is.null(codes)) return(list(d=NA_real_, n_S=0L, n_I=0L, mean_S=NA, mean_I=NA, sd_pooled=NA))

  # Get phenotypes for this trait
  bio_trait <- bio_all[Trait == trait]
  if (nrow(bio_trait) == 0) return(list(d=NA_real_, n_S=0L, n_I=0L, mean_S=NA, mean_I=NA, sd_pooled=NA))

  # Merge SIR codes with phenotypes
  sir_dt <- data.table(Accession=names(codes), SIR=as.character(codes))
  merged <- merge(sir_dt, bio_trait, by="Accession")

  s_accs <- merged[SIR == "S"]
  i_accs <- merged[SIR == "I"]

  if (nrow(s_accs) < 5 || nrow(i_accs) < 5)
    return(list(d=NA_real_, n_S=nrow(s_accs), n_I=nrow(i_accs), mean_S=NA, mean_I=NA, sd_pooled=NA))

  mean_S <- mean(s_accs$Phenotype)
  mean_I <- mean(i_accs$Phenotype)
  sd_S   <- sd(s_accs$Phenotype)
  sd_I   <- sd(i_accs$Phenotype)
  n_S    <- nrow(s_accs)
  n_I    <- nrow(i_accs)

  # Pooled SD
  sd_pooled <- sqrt(((n_S-1)*sd_S^2 + (n_I-1)*sd_I^2) / (n_S+n_I-2))

  # Cohen's d (signed: positive = S superior to I)
  d <- (mean_S - mean_I) / sd_pooled

  list(d=d, n_S=n_S, n_I=n_I, mean_S=mean_S, mean_I=mean_I, sd_pooled=sd_pooled)
}

d_results <- rbindlist(lapply(1:nrow(qtl_trait_pairs), function(i) {
  qtl <- qtl_trait_pairs$QTL[i]; trait <- qtl_trait_pairs$Trait[i]
  r <- compute_d(qtl, trait)
  data.table(QTL=qtl, Trait=trait, d=round(r$d, 3), d_abs=round(abs(r$d), 3),
             n_S=r$n_S, n_I=r$n_I, mean_S=round(r$mean_S, 2), mean_I=round(r$mean_I, 2),
             sd_pooled=round(r$sd_pooled, 3))
}))

cat(sprintf("  Computed: %d with valid d, %d with NA\n",
  sum(!is.na(d_results$d)), sum(is.na(d_results$d))))
cat("\n--- Effect sizes (|d|) by QTL-trait ---\n")
print(d_results[!is.na(d)][order(-d_abs), .(QTL, Trait, d, d_abs, n_S, n_I)], nrows=50)

# =====================================================================
# Step 4: Assign weights to drag pairs
# =====================================================================
cat("\nAssigning pair-level weights...\n")

assign_weight <- function(qtl1, t1, qtl2, t2) {
  d1 <- d_results[QTL==qtl1 & Trait==t1]
  d2 <- d_results[QTL==qtl2 & Trait==t2]
  if (nrow(d1)==0 || nrow(d2)==0) return(NA_real_)
  if (is.na(d1$d) || is.na(d2$d)) return(NA_real_)
  d1$d_abs * d2$d_abs
}

# Internal
int_pairs[, Weight := sapply(1:.N, function(i) assign_weight(QTL[i], Trait1[i], QTL[i], Trait2[i]))]
# Between
bet_pairs[, Weight := sapply(1:.N, function(i) assign_weight(QTL1[i], Trait1[i], QTL2[i], Trait2[i]))]

# Normalize weights: mean weight = 1 for internal, separately for between
int_pairs[, Weight_norm := Weight / mean(Weight, na.rm=TRUE)]
bet_pairs[, Weight_norm := Weight / mean(Weight, na.rm=TRUE)]

cat("\n--- Internal pairs with weights ---\n")
print(int_pairs[, .(QTL, Trait1, Trait2, PureDragPct, Weight=d_results$d_abs[match(paste0(QTL,Trait1),paste0(d_results$QTL,d_results$Trait))] * d_results$d_abs[match(paste0(QTL,Trait2),paste0(d_results$QTL,d_results$Trait))], Weight_norm)][order(-Weight_norm)], nrows=20)

cat("\n--- Between pairs with weights ---\n")
print(bet_pairs[, .(QTL1, QTL2, Trait1, Trait2, PureDragPct, Weight, Weight_norm)][order(-Weight_norm)], nrows=20)

# Save weights
fwrite(int_pairs[, .(QTL, Trait1, Trait2, PureDragPct, Weight, Weight_norm)],
       file.path(OUT, "QTL_internal_pair_weights.csv"), bom=TRUE)
fwrite(bet_pairs[, .(QTL1, QTL2, Trait1, Trait2, PureDragPct, Weight, Weight_norm)],
       file.path(OUT, "QTL_between_pair_weights.csv"), bom=TRUE)
fwrite(d_results, file.path(OUT, "QTL_trait_effect_sizes.csv"), bom=TRUE)

# =====================================================================
# Step 5: Load R-type + meta
# =====================================================================
r_qtl <- fread(file.path(DIR_R, "antagonism_by_QTL.csv"))
r_qtl[, R_Group := fcase(Mean_R >= 20, CAT_HIGH, Mean_R > 0 & Mean_R < 20, CAT_LOW, default = CAT_NONE)]
qtl_r_group <- setNames(r_qtl$R_Group, r_qtl$QTL)

meta <- fread(file.path(META, "yangbenALL.csv"))
meta[, Accession := as.character(Accession)]
meta <- meta[, .(Accession, Type, Date)]

groups <- c("Semi-wild","<1950","1950-1959","1960-1969","1970-1979",
            "1980-1989","1990-1999","2000-2009","2010-2017")

# =====================================================================
# Step 6: Per-accession analysis - BOTH unweighted and weighted
# =====================================================================
cat("\nRunning per-accession drag counting (unweighted + weighted)...\n")

init <- function() list(R_high=0L, R_low=0L,
  int_pure=0L, int_weighted=0, int_valid=0L, int_wsum=0,
  bet_pure=0L, bet_weighted=0, bet_valid=0L, bet_wsum=0)
acc <- setNames(lapply(all_accs, function(a) init()), all_accs)

# R count
for (qtl in names(sir_hash)) {
  grp <- qtl_r_group[qtl]; if (is.na(grp)) grp <- CAT_NONE
  for (trait in names(sir_hash[[qtl]])) {
    codes <- sir_hash[[qtl]][[trait]]
    r_accs <- names(codes)[codes == "R"]
    if (grp == CAT_HIGH) {
      for (a in r_accs) acc[[a]]$R_high <- acc[[a]]$R_high + 1L
    } else if (grp == CAT_LOW) {
      for (a in r_accs) acc[[a]]$R_low <- acc[[a]]$R_low + 1L
    }
  }
}

# Internal drag
for (i in seq_len(nrow(int_pairs))) {
  qtl <- int_pairs$QTL[i]; t1 <- int_pairs$Trait1[i]; t2 <- int_pairs$Trait2[i]
  w  <- int_pairs$Weight_norm[i]
  s1 <- sir_hash[[qtl]][[t1]]; s2 <- sir_hash[[qtl]][[t2]]
  if (is.null(s1) || is.null(s2) || is.na(w)) next
  for (a in all_accs) {
    c1 <- s1[a]; c2 <- s2[a]
    if (!is.na(c1) && !is.na(c2) && c1 != "NA" && c2 != "NA") {
      acc[[a]]$int_valid <- acc[[a]]$int_valid + 1L
      acc[[a]]$int_wsum  <- acc[[a]]$int_wsum + w
      combo <- paste0(c1, "-", c2)
      if (combo %in% PURE_DRAG) {
        acc[[a]]$int_pure     <- acc[[a]]$int_pure + 1L
        acc[[a]]$int_weighted <- acc[[a]]$int_weighted + w
      }
    }
  }
}

# Between drag
for (i in seq_len(nrow(bet_pairs))) {
  q1 <- bet_pairs$QTL1[i]; q2 <- bet_pairs$QTL2[i]
  t1 <- bet_pairs$Trait1[i]; t2 <- bet_pairs$Trait2[i]
  w  <- bet_pairs$Weight_norm[i]
  s1 <- sir_hash[[q1]][[t1]]; s2 <- sir_hash[[q2]][[t2]]
  if (is.null(s1) || is.null(s2) || is.na(w)) next
  for (a in all_accs) {
    c1 <- s1[a]; c2 <- s2[a]
    if (!is.na(c1) && !is.na(c2) && c1 != "NA" && c2 != "NA") {
      acc[[a]]$bet_valid <- acc[[a]]$bet_valid + 1L
      acc[[a]]$bet_wsum  <- acc[[a]]$bet_wsum + w
      combo <- paste0(c1, "-", c2)
      if (combo %in% PURE_DRAG) {
        acc[[a]]$bet_pure     <- acc[[a]]$bet_pure + 1L
        acc[[a]]$bet_weighted <- acc[[a]]$bet_weighted + w
      }
    }
  }
}

acc_df <- rbindlist(lapply(names(acc), function(a) {
  v <- acc[[a]]
  data.table(Accession=a, R_high=v$R_high, R_low=v$R_low, R_total=v$R_high+v$R_low,
    int_pure=v$int_pure, int_wt=v$int_weighted, int_valid=v$int_valid, int_wsum=v$int_wsum,
    bet_pure=v$bet_pure, bet_wt=v$bet_weighted, bet_valid=v$bet_valid, bet_wsum=v$bet_wsum)
}))

# Unweighted rate = drag_pairs / valid_pairs
acc_df[, int_pure_pct := ifelse(int_valid>0, int_pure/int_valid*100, NA)]
acc_df[, bet_pure_pct := ifelse(bet_valid>0, bet_pure/bet_valid*100, NA)]
# Weighted rate = weighted_drag / sum_of_weights
acc_df[, int_wt_pct   := ifelse(int_wsum>0, int_wt/int_wsum*100, NA)]
acc_df[, bet_wt_pct   := ifelse(bet_wsum>0, bet_wt/bet_wsum*100, NA)]

cat(sprintf("\nUnweighted:  int_pure_pct mean=%.1f%%  bet_pure_pct mean=%.1f%%\n",
  mean(acc_df$int_pure_pct, na.rm=TRUE), mean(acc_df$bet_pure_pct, na.rm=TRUE)))
cat(sprintf("Weighted:    int_wt_pct   mean=%.1f%%  bet_wt_pct   mean=%.1f%%\n",
  mean(acc_df$int_wt_pct, na.rm=TRUE), mean(acc_df$bet_wt_pct, na.rm=TRUE)))

# =====================================================================
# Step 7: Merge meta + era
# =====================================================================
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
acc_plot <- acc_all[Group %in% groups]
acc_plot[, Group := factor(Group, levels=groups)]

# =====================================================================
# Step 8: Comparison summary
# =====================================================================
cat("\n\n========== COMPARISON: Unweighted vs Weighted ==========\n")
cat(sprintf("%-14s %12s %12s %8s %12s %12s %8s\n",
  "Era", "int_unwt%", "int_wt%", "d_int", "bet_unwt%", "bet_wt%", "d_bet"))
cat(strrep("-", 85), "\n")
for(g in groups) {
  s <- acc_plot[Group == g]
  iu <- mean(s$int_pure_pct, na.rm=TRUE)
  iw <- mean(s$int_wt_pct, na.rm=TRUE)
  bu <- mean(s$bet_pure_pct, na.rm=TRUE)
  bw <- mean(s$bet_wt_pct, na.rm=TRUE)
  cat(sprintf("%-14s %11.1f%% %11.1f%% %+7.1f  %11.1f%% %11.1f%% %+7.1f\n",
    g, iu, iw, iw-iu, bu, bw, bw-bu))
}

# Absolute count comparison
cat("\n\n========== COMPARISON: Absolute Counts (Unweighted vs Weighted) ==========\n")
cat(sprintf("%-14s %12s %12s %8s %12s %12s %8s\n",
  "Era", "int_count", "int_wt", "d_int", "bet_count", "bet_wt", "d_bet"))
cat(strrep("-", 85), "\n")
for(g in groups) {
  s <- acc_plot[Group == g]
  iu <- mean(s$int_pure)
  iw <- mean(s$int_wt)
  bu <- mean(s$bet_pure)
  bw <- mean(s$bet_wt)
  cat(sprintf("%-14s %11.2f %11.2f %+7.2f  %11.2f %11.2f %+7.2f\n",
    g, iu, iw, iw-iu, bu, bw, bw-bu))
}

# =====================================================================
# Step 9: KW + Dunn's on weighted rates
# =====================================================================
fmt_p <- function(p) {
  if (p < 0.0001) return("p < 0.0001")
  if (p < 0.001)  return(sprintf("p = %.4f", p))
  if (p < 0.01)   return(sprintf("p = %.4f", p))
  if (p < 0.05)   return(sprintf("p = %.4f", p))
  return(sprintf("p = %.3f (ns)", p))
}

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

cat("\n========== KW: Weighted drag rates ==========\n")
for (yv in c("int_wt_pct","bet_wt_pct")) {
  kw <- kruskal.test(as.formula(paste0(yv, " ~ Group")), data=acc_plot)
  cat(sprintf("%s: %s\n", yv, fmt_p(kw$p.value)))
}

# =====================================================================
# Step 10: Plots - weighted vs unweighted
# =====================================================================
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

make_cmp_plot <- function(y_unwt, y_wt, title_text, y_label) {
  d1 <- acc_plot[, .(Group, Rate=get(y_unwt))]
  d1[, Method := "Unweighted"]
  d2 <- acc_plot[, .(Group, Rate=get(y_wt))]
  d2[, Method := "Weighted (|d1 x d2|)"]
  dm <- rbind(d1, d2)
  dm[, Method := factor(Method, levels=c("Unweighted","Weighted (|d1 x d2|)"))]

  # Compute means + letters per method
  m_list <- list()
  for (mth in unique(dm$Method)) {
    sub <- dm[Method == mth]
    sub[, Group := factor(Group, levels=groups)]
    kw <- kruskal.test(Rate ~ Group, data=sub)
    mm <- sub[, .(Mean=mean(Rate, na.rm=TRUE), yMax=max(Rate, na.rm=TRUE)), by=Group]
    ll <- get_letters(sub, "Rate")
    mm <- merge(mm, ll, by="Group")
    yrange <- diff(range(c(0, mm$yMax)))
    if (yrange == 0 || is.na(yrange)) yrange <- 1
    mm[, LabelY := yMax + yrange * 0.06]
    mm[, Method := mth]
    mm[, KW_p := kw$p.value]
    m_list[[mth]] <- mm
  }
  m_all <- rbindlist(m_list)

  ymax <- max(dm$Rate, na.rm=TRUE)
  ggplot(dm, aes(x=Group, y=Rate, fill=Group)) +
    geom_boxplot(outlier.size=0.3, alpha=0.85, width=0.6) +
    geom_label(data=m_all, aes(x=Group, y=LabelY, label=sprintf("%.1f %s", Mean, Letter)),
               fill="white", alpha=0.85, size=2.3, fontface="bold",
               vjust=-0.5, label.padding=unit(0.1, "lines")) +
    facet_wrap(~Method, ncol=1, scales="free_y") +
    scale_fill_manual(values=group_cols) + scale_x_discrete(labels=n_vec) +
    labs(x="", y=y_label, title=title_text) + white_theme
}

p_int <- make_cmp_plot("int_pure_pct", "int_wt_pct",
  "QTL-internal drag: Unweighted vs Phenotype-weighted", "Internal drag rate (%)")
p_bet <- make_cmp_plot("bet_pure_pct", "bet_wt_pct",
  "QTL-between drag: Unweighted vs Phenotype-weighted", "Between drag rate (%)")

ggsave(file.path(OUT, "internal_weighted_vs_unweighted.pdf"), p_int, width=10, height=7.5)
ggsave(file.path(OUT, "internal_weighted_vs_unweighted.png"), p_int, width=10, height=7.5, dpi=300)
ggsave(file.path(OUT, "between_weighted_vs_unweighted.pdf"), p_bet, width=10, height=7.5)
ggsave(file.path(OUT, "between_weighted_vs_unweighted.png"), p_bet, width=10, height=7.5, dpi=300)

# Absolute count plots
p_int_abs <- make_cmp_plot("int_pure", "int_wt",
  "QTL-internal drag (absolute count): Unweighted vs Weighted", "Internal drag count")
p_bet_abs <- make_cmp_plot("bet_pure", "bet_wt",
  "QTL-between drag (absolute count): Unweighted vs Weighted", "Between drag count")

ggsave(file.path(OUT, "internal_absolute_count.pdf"), p_int_abs, width=10, height=7.5)
ggsave(file.path(OUT, "internal_absolute_count.png"), p_int_abs, width=10, height=7.5, dpi=300)
ggsave(file.path(OUT, "between_absolute_count.pdf"), p_bet_abs, width=10, height=7.5)
ggsave(file.path(OUT, "between_absolute_count.png"), p_bet_abs, width=10, height=7.5, dpi=300)

# =====================================================================
# Step 11: Save
# =====================================================================
summ_rate <- acc_plot[, .(
  N = .N,
  int_unwt_pct = round(mean(int_pure_pct, na.rm=TRUE), 1),
  int_wt_pct   = round(mean(int_wt_pct, na.rm=TRUE), 1),
  bet_unwt_pct = round(mean(bet_pure_pct, na.rm=TRUE), 1),
  bet_wt_pct   = round(mean(bet_wt_pct, na.rm=TRUE), 1)
), by=Group][order(Group)]

summ_abs <- acc_plot[, .(
  N = .N,
  int_unwt = round(mean(int_pure), 2),
  int_wt   = round(mean(int_wt), 2),
  bet_unwt = round(mean(bet_pure), 2),
  bet_wt   = round(mean(bet_wt), 2)
), by=Group][order(Group)]

cat("\n\n========== FINAL SUMMARY (Rates) ==========\n")
print(summ_rate)
cat("\n\n========== FINAL SUMMARY (Absolute counts) ==========\n")
print(summ_abs)

fwrite(summ_rate, file.path(OUT, "weighted_vs_unweighted_by_period_rates.csv"), bom=TRUE)
fwrite(summ_abs, file.path(OUT, "weighted_vs_unweighted_by_period_counts.csv"), bom=TRUE)
fwrite(acc_all, file.path(OUT, "three_types_per_accession_weighted.csv"), bom=TRUE)

cat("\n=== Done ===\n")
cat("Output:", OUT, "\n")
