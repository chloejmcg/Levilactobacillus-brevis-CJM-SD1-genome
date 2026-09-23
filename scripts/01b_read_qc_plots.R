#!/usr/bin/env Rscript
# 01b_read_qc_plots.R — summarise and plot per-read length/quality
# Usage: Rscript 01b_read_qc_plots.R <per_read_len_qual.tsv> <outdir>
suppressPackageStartupMessages(library(ggplot2))

args   <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outdir <- args[2]
genome_size <- 2.5e6   # from the Plasmidsaurus assembly

d <- read.delim(infile, header = FALSE, quote = "")[, 1:3]
names(d) <- c("read_id", "length", "mean_q")

# N50: the length L such that reads >= L contain half of all bases
n50 <- function(x) {
  x <- sort(x, decreasing = TRUE)
  x[which(cumsum(as.numeric(x)) >= sum(as.numeric(x)) / 2)[1]]
}

total_bp <- sum(as.numeric(d$length))
summ <- data.frame(
  metric = c("reads", "total_bases", "mean_length", "read_N50", "max_length",
             "median_mean_Q", "reads_Q10_plus", "reads_Q20_plus", "est_coverage_x"),
  value  = c(nrow(d), total_bp, round(mean(d$length)), n50(d$length), max(d$length),
             round(median(d$mean_q), 1), sum(d$mean_q >= 10), sum(d$mean_q >= 20),
             round(total_bp / genome_size, 1))
)
write.table(summ, file.path(outdir, "read_qc_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
print(summ, row.names = FALSE)

p1 <- ggplot(d, aes(length, mean_q)) +
  geom_bin2d(bins = 100) +
  scale_x_log10() +
  scale_fill_viridis_c(trans = "log10") +
  labs(x = "Read length (bp, log scale)", y = "Mean read quality (Phred)",
       fill = "Reads", title = "DLVP55_2 raw ONT reads") +
  theme_bw()
ggsave(file.path(outdir, "length_vs_quality.png"), p1, width = 6, height = 4.5, dpi = 300)

p2 <- ggplot(d, aes(length)) +
  geom_histogram(bins = 100, fill = "grey40") +
  geom_vline(xintercept = n50(d$length), linetype = "dashed", colour = "red") +
  scale_x_log10() +
  labs(x = "Read length (bp, log scale)", y = "Reads",
       title = paste0("Read length distribution (N50 = ", n50(d$length), " bp)")) +
  theme_bw()
ggsave(file.path(outdir, "length_histogram.png"), p2, width = 6, height = 4.5, dpi = 300)
