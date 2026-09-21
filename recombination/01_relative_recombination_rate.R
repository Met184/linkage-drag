#!/usr/bin/env Rscript
# 01_relative_recombination_rate.R
# Per-QTL absolute and relative recombination rate in the three populations.
# Usage: Rscript 01_relative_recombination_rate.R [--work DIR] [--out DIR]

library(data.table)

# ---- CONFIG -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1L) args[i + 1L] else default
}

WORK     <- get_arg("--work", ".")
MAP_DIR  <- get_arg("--maps", file.path(WORK, "maps"))
QTL_FILE <- get_arg("--qtl",  file.path(WORK, "QTL_intervals.csv"))
DRAG_DIR <- get_arg("--drag", file.path(WORK, "drag"))
OUT_DIR  <- get_arg("--out",  file.path(WORK, "relative_rate"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("Recombination rate of three populations - raw values + relative rate\n")
cat("============================================================\n\n")

# ====== Step 1: build the raw recombination maps ======
cat("Step 1: building the raw recombination maps\n")

# Group-D (already in the expected units)
gD <- fread(file.path(MAP_DIR, "Group_D_published_cM_per_Mb.txt"), header = FALSE)
setnames(gD, c("Chr", "Start", "End", "Recomb_Rate"))
gD[, `:=`(Chr = as.character(Chr), Start = as.numeric(Start),
           End = as.numeric(End), Recomb_Rate = as.numeric(Recomb_Rate))]

# Group-C raw (Rho / 50)
gC_raw <- fread(file.path(MAP_DIR, "Group_C_fastEPRR_rho_50kb.txt"))
setnames(gC_raw, c("Chr", "Start", "End", "Rho", "CIL", "CIR"))
gC_raw[, `:=`(Chr = as.character(gsub("^chr", "", Chr)),
               Start = as.numeric(Start), End = as.numeric(End), Rho = as.numeric(Rho))]
gC <- gC_raw[, .(Chr, Start, End, Recomb_Rate = Rho / 50)]
gC <- gC[!is.na(Recomb_Rate) & !is.na(Start) & !is.na(End)]

# Group-B raw (Rho / 50)
gB_raw <- fread(file.path(MAP_DIR, "Group_B_fastEPRR_rho_50kb.txt"))
setnames(gB_raw, c("Chr", "Start", "End", "Rho", "CIL", "CIR"))
gB_raw[, `:=`(Chr = as.character(gsub("^chr", "", Chr)),
               Start = as.numeric(Start), End = as.numeric(End), Rho = as.numeric(Rho))]
gB <- gB_raw[, .(Chr, Start, End, Recomb_Rate = Rho / 50)]
gB <- gB[!is.na(Recomb_Rate) & !is.na(Start) & !is.na(End)]

maps <- list(Group_D = gD, Group_C = gC, Group_B = gB)

# Genome-wide statistics
gstat <- rbindlist(lapply(names(maps), function(p) {
  x <- maps[[p]]$Recomb_Rate
  data.table(Population = p, N_windows = length(x), Mean = mean(x),
             Median = median(x), SD = sd(x), CV = sd(x)/mean(x),
             Min = min(x), Max = max(x))
}))
cat("Genome-wide statistics:\n")
print(gstat)

# ====== Step 2: recombination rate of each QTL ======
cat("\nStep 2: recombination rate of each QTL\n")

qtl_pos <- fread(QTL_FILE)

calc_rec <- function(chr, s, e, map) {
  cm <- map[Chr == chr]
  if (nrow(cm) == 0) return(NA_real_)
  ov <- cm[Start <= e & End >= s]
  if (nrow(ov) == 0) return(NA_real_)
  tw <- 0; tv <- 0
  for (j in 1:nrow(ov)) {
    ol <- min(e, ov$End[j]) - max(s, ov$Start[j])
    if (ol > 0) { tw <- tw + ol; tv <- tv + ov$Recomb_Rate[j] * ol }
  }
  if (tw > 0) tv/tw else NA_real_
}

qr <- data.table()
for (i in 1:nrow(qtl_pos)) {
  chr <- qtl_pos$chr[i]; s <- qtl_pos$new_QTL_start[i]; e <- qtl_pos$new_QTL_end[i]
  qr <- rbind(qr, data.table(
    QTL_id = qtl_pos$QTL_id[i], Chromosome = chr,
    Rec_Group_D = calc_rec(chr, s, e, gD),
    Rec_Group_C = calc_rec(chr, s, e, gC),
    Rec_Group_B = calc_rec(chr, s, e, gB)
  ))
}
# Relative recombination rate
qr[, `:=`(Rel_Group_D = Rec_Group_D / gstat[Population == "Group_D", Mean],
          Rel_Group_C = Rec_Group_C / gstat[Population == "Group_C", Mean],
          Rel_Group_B = Rec_Group_B / gstat[Population == "Group_B", Mean])]
cat(sprintf("computed %d QTLs\n", nrow(qr)))

# ====== Step 3: QTL classification ======
cat("\nStep 3: QTL classification\n")

antago    <- fread(file.path(DRAG_DIR, "antagonism_QTL_trait.csv"))
intra_drag <- fread(file.path(DRAG_DIR, "QTL_internal_drag.csv"))
inter_drag <- fread(file.path(DRAG_DIR, "QTL_between_drag.csv"))
# The input table must use the English column names:
#   Trait1 Trait2 Type N SS II SI IS R_total Synergy PureDrag RDrag
#   SynergyPct PureDragPct RDragPct Category

# Drag QTL definition: R > 10% or PureDrag > 10% in any analysis
drag_qtls <- unique(unlist(c(
  antago[RPct > 10, QTL],
  intra_drag[PureDragPct > 10, QTL],
  inter_drag[PureDragPct > 10, c(QTL1, QTL2)]
)))
drag_qtls <- intersect(drag_qtls, qr[!is.na(Rec_Group_D), QTL_id])
nondrag_qtls <- setdiff(qr[!is.na(Rec_Group_D), QTL_id], drag_qtls)
cat(sprintf("Drag: %d, Non-drag: %d, Total: %d\n", length(drag_qtls), length(nondrag_qtls),
            length(drag_qtls) + length(nondrag_qtls)))

qr[, Group := ifelse(QTL_id %in% drag_qtls, "Drag", "Non-drag")]

# ====== Step 4: group statistics ======
cat("\nStep 4: group statistics\n")

gs <- data.table()
for (pop in c("Group_D", "Group_C", "Group_B")) {
  col_abs <- paste0("Rec_", pop); col_rel <- paste0("Rel_", pop)
  gm <- gstat[Population == pop, Mean]

  for (grp in c("Drag", "Non-drag")) {
    vals_abs <- qr[Group == grp, get(col_abs)]; vals_abs <- vals_abs[!is.na(vals_abs)]
    vals_rel <- qr[Group == grp, get(col_rel)]; vals_rel <- vals_rel[!is.na(vals_rel)]
    gs <- rbind(gs, data.table(Population = pop, Group = grp,
      N = length(vals_abs), Mean_abs = mean(vals_abs), SD_abs = sd(vals_abs),
      Mean_rel = mean(vals_rel), SD_rel = sd(vals_rel)))
  }
  # statistical tests
  d_abs <- qr[Group == "Drag", get(col_abs)]; d_abs <- d_abs[!is.na(d_abs)]
  n_abs <- qr[Group == "Non-drag", get(col_abs)]; n_abs <- n_abs[!is.na(n_abs)]
  cat(sprintf("  %s: Drag=%.2f, Non-drag=%.2f | DvsND p=%.4f | DvsGenome p=%.4f\n",
    pop, mean(d_abs), mean(n_abs),
    wilcox.test(d_abs, n_abs)$p.value, wilcox.test(d_abs, mu = gm)$p.value))
}

cat("\nGroup statistics table:\n")
print(gs)

# ====== Step 5: antagonism R > 10% ======
cat("\nStep 5: antagonism R > 10%\n")

calc_pheno <- function(qtl_id, trait) {
  qnum <- gsub("QTL", "", qtl_id)
  ppf <- file.path(DRAG_DIR, "per_pair", sprintf("QTL_%s_%s_SIR.csv", qnum, trait))
  if (!file.exists(ppf)) return(list(n_S=NA, n_I=NA, d=NA, mean_S=NA, mean_I=NA))
  pp <- fread(ppf); setnames(pp, c("Accession", "SIR"))
  bf <- file.path(DRAG_DIR, "phenotypes", paste0(trait, ".csv"))
  if (!file.exists(bf)) return(list(n_S=NA, n_I=NA, d=NA, mean_S=NA, mean_I=NA))
  bio <- fread(bf); setnames(bio, c("SRR", "bg", "pheno"))
  bio[, pheno := as.numeric(pheno)]; bio <- bio[!is.na(pheno)]
  m <- merge(pp, bio, by.x = "Accession", by.y = "SRR")
  sv <- m[SIR == "S", pheno]; iv <- m[SIR == "I", pheno]
  ns <- length(sv); ni <- length(iv)
  if (ns < 2 || ni < 2) return(list(n_S=ns, n_I=ni, d=NA, mean_S=if(ns>0) mean(sv) else NA, mean_I=if(ni>0) mean(iv) else NA))
  ps <- sqrt(((ns-1)*var(sv) + (ni-1)*var(iv)) / (ns+ni-2))
  list(n_S=ns, n_I=ni, d=(mean(sv)-mean(iv))/ps, mean_S=mean(sv), mean_I=mean(iv))
}

ah <- antago[RPct > 10][order(-RPct)]
ar <- data.table()
for (i in 1:nrow(ah)) {
  r <- ah[i]; qtl <- r$QTL; trait <- r$Trait
  rec <- qr[QTL_id == qtl]
  ph <- calc_pheno(qtl, trait)
  ar <- rbind(ar, data.table(
    QTL = qtl, Trait = trait, R_pct = r$RPct, S_pct = r$SPct, I_pct = r$IPct,
    Rec_Group_D = rec$Rec_Group_D, Rel_Group_D = rec$Rel_Group_D,
    Rec_Group_C = rec$Rec_Group_C, Rel_Group_C = rec$Rel_Group_C,
    Rec_Group_B = rec$Rec_Group_B, Rel_Group_B = rec$Rel_Group_B,
    n_S = ph$n_S, mean_S = ph$mean_S, n_I = ph$n_I, mean_I = ph$mean_I, cohens_d = ph$d
  ), fill = TRUE)
}
cat(sprintf("R>10%%: %d cases, with d: %d\n", nrow(ar), sum(!is.na(ar$cohens_d))))
print(ar[, .(QTL, Trait, R_pct, Rec_Group_D, Rel_Group_D, Rec_Group_C, Rel_Group_C, cohens_d)], nrows=99)

# ====== Step 6: QTL-internal drag ======
cat("\nStep 6: QTL-internal drag, PureDrag > 10%\n")

ih <- intra_drag[PureDragPct > 10][order(-PureDragPct)]
br1 <- data.table()
for (i in 1:nrow(ih)) {
  r <- ih[i]; rec <- qr[QTL_id == r$QTL]
  if (nrow(rec) == 0) next
  br1 <- rbind(br1, data.table(
    QTL = r$QTL, Pair = paste0(r$Trait1, "x", r$Trait2),
    N = r$N, PureDragPct = r$PureDragPct, Category = r$Category,
    Rec_Group_D = rec$Rec_Group_D, Rel_Group_D = rec$Rel_Group_D,
    Rec_Group_C = rec$Rec_Group_C, Rel_Group_C = rec$Rel_Group_C,
    Rec_Group_B = rec$Rec_Group_B, Rel_Group_B = rec$Rel_Group_B
  ))
}
gD_mean <- gstat[Population == "Group_D", Mean]
cat(sprintf("%d cases in total, all below the Group_D mean (%.2f): %s\n",
  nrow(br1), gD_mean,
  ifelse(all(br1$Rec_Group_D < gD_mean, na.rm = TRUE), "YES", "NO")))
print(br1[, .(QTL, Pair, PureDragPct, Rec_Group_D, Rel_Group_D, Rec_Group_C, Rel_Group_C)], nrows=99)

# ====== Step 7: QTL-between drag ======
cat("\nStep 7: QTL-between drag, PureDrag > 10%\n")

ieh <- inter_drag[PureDragPct > 10][order(-PureDragPct)]
br2 <- data.table()
for (i in 1:nrow(ieh)) {
  r <- ieh[i]; q1 <- r$QTL1; q2 <- r$QTL2
  p1 <- qtl_pos[QTL_id == q1]; p2 <- qtl_pos[QTL_id == q2]
  if (nrow(p1) == 0 || nrow(p2) == 0) next

  chr <- p1$chr[1]
  if (chr != p2$chr[1]) next

  ws <- min(p1$new_QTL_start, p2$new_QTL_start)
  we <- max(p1$new_QTL_end, p2$new_QTL_end)
  if (p1$new_QTL_end < p2$new_QTL_start) { bs <- p1$new_QTL_end; be <- p2$new_QTL_start }
  else { bs <- p2$new_QTL_end; be <- p1$new_QTL_start }

  bgD <- calc_rec(chr, bs, be, gD); bgC <- calc_rec(chr, bs, be, gC); bgB <- calc_rec(chr, bs, be, gB)

  br2 <- rbind(br2, data.table(
    QTL1 = q1, QTL2 = q2, Pair = paste0(r$Trait1, "x", r$Trait2),
    N = r$N, PureDragPct = r$PureDragPct, PhysDist_Mb = abs(be-bs)/1e6,
    Between_Group_D = bgD, Rel_Between_Group_D = bgD / gstat[Population == "Group_D", Mean],
    Between_Group_C = bgC, Rel_Between_Group_C = bgC / gstat[Population == "Group_C", Mean],
    Between_Group_B = bgB, Rel_Between_Group_B = bgB / gstat[Population == "Group_B", Mean]
  ))
}
print(br2[, .(QTL1, QTL2, Pair, PureDragPct, Between_Group_D, Rel_Between_Group_D,
              Between_Group_C, Rel_Between_Group_C)], nrows=99)

# ====== Step 8: Ghir_A07 tandem array ======
cat("\nStep 8: Ghir_A07 QTL051-052-053-054\n")

cat("Within each QTL:\n")
ca07 <- qr[QTL_id %in% c("QTL051","QTL052","QTL053","QTL054")]
ca07 <- merge(ca07, qtl_pos[, .(QTL_id, new_QTL_start, new_QTL_end)], by = "QTL_id")
ca07[, `:=`(Len_kb = round((new_QTL_end - new_QTL_start)/1e3, 1))]
print(ca07[, .(QTL_id, Len_kb, Rec_Group_D, Rel_Group_D, Rec_Group_C, Rel_Group_C, Rec_Group_B, Rel_Group_B)])

cat("Between pairs:\n")
cqtls <- c("QTL051","QTL052","QTL053","QTL054"); cp <- t(combn(cqtls, 2))
cb <- data.table()
for (k in 1:nrow(cp)) {
  q1 <- cp[k,1]; q2 <- cp[k,2]
  p1 <- qtl_pos[QTL_id == q1]; p2 <- qtl_pos[QTL_id == q2]
  chr <- p1$chr[1]
  if (p1$new_QTL_end < p2$new_QTL_start) { bs <- p1$new_QTL_end; be <- p2$new_QTL_start }
  else { bs <- p2$new_QTL_end; be <- p1$new_QTL_start }
  cb <- rbind(cb, data.table(QTL1 = q1, QTL2 = q2, Dist_Mb = round(abs(be-bs)/1e6, 3),
    B_Group_D = calc_rec(chr, bs, be, gD), B_Group_D_rel = calc_rec(chr, bs, be, gD) / gstat[Population == "Group_D", Mean],
    B_Group_C = calc_rec(chr, bs, be, gC), B_Group_C_rel = calc_rec(chr, bs, be, gC) / gstat[Population == "Group_C", Mean],
    B_Group_B = calc_rec(chr, bs, be, gB), B_Group_B_rel = calc_rec(chr, bs, be, gB) / gstat[Population == "Group_B", Mean]))
}
print(cb)

# ====== Step 9: save results ======
cat("\nStep 9: saving\n")

fwrite(gstat, file.path(OUT_DIR, "01_genomewide_stats.csv"))
fwrite(gs,    file.path(OUT_DIR, "02_drag_vs_nondrag.csv"))
fwrite(ar,    file.path(OUT_DIR, "03_antagonism_R_gt10.csv"))
fwrite(br1,   file.path(OUT_DIR, "04_QTL_internal_drag.csv"))
fwrite(br2,   file.path(OUT_DIR, "05_QTL_between_drag.csv"))
fwrite(ca07,  file.path(OUT_DIR, "06_GhirA07_QTLs.csv"))
fwrite(cb,    file.path(OUT_DIR, "06_GhirA07_QTL_pairs.csv"))
fwrite(qr,    file.path(OUT_DIR, "07_all_QTL_rates.csv"))

# ====== summary report ======
sink(file.path(OUT_DIR, "analysis_report.txt"))
cat("============================================================\n")
cat("Recombination rate of three populations - raw values (Rho/50, uncalibrated)\n")
cat("Cross-population comparison: relative rate = QTL value / genome-wide mean\n")
cat("============================================================\n\n")

cat("[Genome-wide statistics]\n")
print(gstat)

cat("\n[Group statistics]\n")
cat("Mean_rel = 1.0 equals the genome-wide mean; < 1.0 is below the mean (cold region)\n\n")
print(gs)

cat("\n[Statistical tests]\n")
for (pop in c("Group_D", "Group_C", "Group_B")) {
  ca <- paste0("Rec_", pop)
  da <- qr[Group == "Drag", get(ca)]; da <- da[!is.na(da)]
  na <- qr[Group == "Non-drag", get(ca)]; na <- na[!is.na(na)]
  gm <- gstat[Population == pop, Mean]
  cat(sprintf("%s: Drag %.2f vs Non-drag %.2f, p=%.4f | Drag vs Genome(%.2f), p=%.4f\n",
    pop, mean(da), mean(na), wilcox.test(da, na)$p.value, gm, wilcox.test(da, mu=gm)$p.value))
}

cat("\n\n[Antagonism R > 10%]\n")
cat(sprintf("%d cases\n", nrow(ar)))
cat(sprintf("Group_D absolute range: %.2f-%.2f  relative range: %.2f-%.2f\n",
  min(ar$Rec_Group_D, na.rm=TRUE), max(ar$Rec_Group_D, na.rm=TRUE),
  min(ar$Rel_Group_D, na.rm=TRUE), max(ar$Rel_Group_D, na.rm=TRUE)))
cat(sprintf("|d| range: %.2f-%.2f\n", min(abs(ar$cohens_d), na.rm=TRUE), max(abs(ar$cohens_d), na.rm=TRUE)))

cat("\n[QTL-internal drag]\n")
cat(sprintf("%d cases, all below the Group_D mean (%.2f): %s\n",
  nrow(br1), gD_mean, ifelse(all(br1$Rec_Group_D < gD_mean, na.rm=TRUE), "yes", "no")))
cat(sprintf("Group_D range: %.2f-%.2f  relative: %.2f-%.2f\n",
  min(br1$Rec_Group_D, na.rm=TRUE), max(br1$Rec_Group_D, na.rm=TRUE),
  min(br1$Rel_Group_D, na.rm=TRUE), max(br1$Rel_Group_D, na.rm=TRUE)))

cat("\n[QTL-between drag]\n")
cat(sprintf("%d QTL pairs\n", nrow(br2)))
cat(sprintf("Between Group_D range: %.2f-%.2f  relative: %.2f-%.2f\n",
  min(br2$Between_Group_D, na.rm=TRUE), max(br2$Between_Group_D, na.rm=TRUE),
  min(br2$Rel_Between_Group_D, na.rm=TRUE), max(br2$Rel_Between_Group_D, na.rm=TRUE)))

cat("\n[Ghir_A07 tandem array]\n")
cat(sprintf("QTL051-054 overall: Group_D=%.2f (%.2fx), Group_C=%.2f (%.2fx), Group_B=%.2f (%.2fx)\n",
  mean(ca07$Rec_Group_D), mean(ca07$Rel_Group_D),
  mean(ca07$Rec_Group_C), mean(ca07$Rel_Group_C),
  mean(ca07$Rec_Group_B), mean(ca07$Rel_Group_B)))

sink()

cat("\n===== done =====\n")
cat("output: ", OUT_DIR, "\n")
