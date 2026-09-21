#!/usr/bin/env Rscript
# 02_cross_population_correlation.R
# Pairwise Spearman correlation of the relative rates on the common 100-kb windows.
# Usage: Rscript 02_cross_population_correlation.R [--work DIR]

library(data.table)

# ---- CONFIG -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1L) args[i + 1L] else default
}

WORK    <- get_arg("--work", ".")
MAP_DIR <- get_arg("--maps", file.path(WORK, "maps"))
OUT_DIR <- get_arg("--out",  WORK)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

gD <- fread(file.path(MAP_DIR, "Group_D_recombination_cM_per_Mb.csv"))
gB <- fread(file.path(MAP_DIR, "Group_B_recombination_cM_per_Mb.csv"))
gC <- fread(file.path(MAP_DIR, "Group_C_recombination_cM_per_Mb.csv"))
setnames(gD, c("Chr","Start","End","Rate"))
setnames(gB, c("Chr","Start","End","Rate"))
setnames(gC, c("Chr","Start","End","Rate"))

BIN <- 100000
all_chrs <- sort(unique(c(gD$Chr, gB$Chr, gC$Chr)))

bins <- data.table()
for (chr in all_chrs) {
  for (pop in c("Group_D","Group_B","Group_C")) {
    map <- switch(pop, Group_D=gD, Group_B=gB, Group_C=gC)
    cm <- map[Chr==chr]; if(nrow(cm)==0) next
    mx <- max(cm$End)
    brks <- seq(0, mx+BIN, by=BIN)
    cm[, bin:=findInterval(Start, brks)]
    bin_rates <- cm[,.(Rate=sum(Rate*(End-Start))/sum(End-Start)), by=bin]
    bin_rates[,`:=`(Chr=chr, Pop=pop)]
    bins <- rbind(bins, bin_rates, fill=TRUE)
  }
}

bw <- dcast(bins, Chr+bin~Pop, value.var="Rate")
bw <- bw[!is.na(Group_D) & !is.na(Group_B) & !is.na(Group_C)]

# genome-wide mean rate, then convert to relative values
genome_mean_gD <- sum(bw$Group_D)     / nrow(bw)
genome_mean_gB <- sum(bw$Group_B) / nrow(bw)
genome_mean_gC <- sum(bw$Group_C) / nrow(bw)

bw[, Rel_Group_D     := Group_D     / genome_mean_gD]
bw[, Rel_Group_B := Group_B / genome_mean_gB]
bw[, Rel_Group_C := Group_C / genome_mean_gC]

cat(sprintf("Genome mean (relative units): Group_D=%.2f, Group_B=%.2f, Group_C=%.2f\n\n",
  genome_mean_gD, genome_mean_gB, genome_mean_gC))

COLORS <- c("white","grey80","blue","darkblue")

plot_pair <- function(x, y, xlab, ylab, fn) {
  rp <- round(cor(x, y, method="pearson"), 3)
  rs <- round(cor(x, y, method="spearman"), 3)

  pdf(file.path(OUT_DIR, paste0(fn, ".pdf")), width=7, height=7)
  par(mar=c(4.5,4.5,2,2))
  smoothScatter(x, y, xlab=xlab, ylab=ylab, main="",
    nbin=512, nrpoints=0,
    colramp=colorRampPalette(COLORS))
  abline(0, 1, lty=2, col="red")
  legend("topleft",
    sprintf("Spearman rho = %.3f", rs),
    bty="n", text.col="black", cex=1.2)
  dev.off()

  png(file.path(OUT_DIR, paste0(fn, ".png")), width=3000, height=3000, res=450)
  par(mar=c(4.5,4.5,2,2))
  smoothScatter(x, y, xlab=xlab, ylab=ylab, main="",
    nbin=512, nrpoints=0,
    colramp=colorRampPalette(COLORS))
  abline(0, 1, lty=2, col="red")
  legend("topleft",
    sprintf("Spearman rho = %.3f", rs),
    bty="n", text.col="black", cex=1.4)
  dev.off()

  cat(sprintf("%s vs %s: Spearman rho=%.3f\n", xlab, ylab, rs))
}

plot_pair(bw$Rel_Group_D,     bw$Rel_Group_B, "Group-D Rel", "Group-B Rel", "Fig_Genomewide_Rel_Correlation_Group_D_vs_Group_B")
plot_pair(bw$Rel_Group_D,     bw$Rel_Group_C, "Group-D Rel", "Group-C Rel", "Fig_Genomewide_Rel_Correlation_Group_D_vs_Group_C")
plot_pair(bw$Rel_Group_B, bw$Rel_Group_C, "Group-B Rel", "Group-C Rel", "Fig_Genomewide_Rel_Correlation_Group_B_vs_Group_C")

cat("\nDone. N =", nrow(bw), "windows\n")
