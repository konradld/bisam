# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     Figure: what the non-local prior buys (rates of evidence accumulation)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Consumes the cache from scripts/16_prior_comparison_rates.R and draws
#   output/simulation/figure_prior_rates.pdf
#
# Include in the paper with:
#   \includegraphics[width=\linewidth]{img/figure_prior_rates.pdf}
# NO scale= factor and no .9 — the figure is drawn at final print size.
#
# Suggested caption (the figure deliberately carries no titles):
#   "Rates of evidence accumulation under local and non-local slab priors, on one
#    unit of the step-indicator design with candidate dates 3..T-1. Panel (a):
#    mean log Bayes factor for including a spurious break when no break is
#    present. Panel (b): mean log Bayes factor for the true break when a break of
#    three error standard deviations is present. Panel (c): evidence against one
#    spurious break, |log BF|, against the multiplicity burden log K implied by
#    K = N(T-3) candidate dates. All priors are calibrated to
#    P(|gamma| <= sigma | tau) = 0.01 as in Section 2.3; the Zellner arm uses the
#    conventional unit-information g = n, and the two local arms coincide in slope
#    despite very different tau, confirming the rate is calibration-invariant."
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
source("./scripts/journal_style.R")

CACHE <- "./output/simulation/prior_rate_comparison.RDS"
OUT   <- "./output/simulation/figure_prior_rates.pdf"

if (!file.exists(CACHE))
  stop("Run scripts/16_prior_comparison_rates.R first to build ", CACHE)

z    <- readRDS(CACHE)
res  <- z$res
Ni   <- z$Ni

# Plot order: locals first (the baseline being displaced), piMOM last so it draws
# on top. Two local families share one colour: they are one family, and the point
# is the family's rate, not the individual prior.
arms <- list(
  list(key = "zellner_g_n",  col = JCOL$local, lty = 1, lab = "Local: Zellner"),
  list(key = "normalid_cal", col = JCOL$local, lty = 3, lab = "Local: normal slab"),
  list(key = "pmom",         col = JCOL$mom,   lty = 2, lab = "pMOM"),
  list(key = "pimom",        col = JCOL$imom,  lty = 1, lab = "piMOM")
)

get <- function(key, arm) {
  d <- res[res$key == key & res$arm == arm, ]
  d[order(d$T), ]
}

TS <- sort(unique(res$T))
xt <- c(20, 80, 320, 640)

jrnl_pdf(OUT, height = 2.55)
jrnl_par(mfrow = c(1, 3), mar = c(2.9, 3.0, 1.2, 1.6))

# ---- (a) evidence AGAINST a spurious break ----------------------------------
# The paper's cut claim (tex:270-276) lives here: locals crawl at -1/2 log T while
# piMOM plunges root-exponentially.
# x is on a log scale throughout, which LINEARISES the polynomial rates: a local
# prior is a straight line of slope -1/2 and pMOM a straight line of slope -3/2,
# while piMOM is not a line at all -- it accelerates away. That qualitative
# difference (line vs. curve) is the whole point of the panel.
h0 <- lapply(arms, function(a) get(a$key, "H0"))
lo_a <- min(unlist(lapply(h0, function(d) d$mean)))
yt_a <- seq(-50, 0, 10)
ylim_a <- c(min(lo_a * 1.06, -10), 2)

plot(NA, xlim = range(TS), ylim = ylim_a, xlab = "Time periods (T)",
     ylab = "log Bayes factor", main = "", axes = FALSE, log = "x")
jrnl_grid(h = yt_a)
abline(h = 0, col = JCOL$mute, lwd = JRNL$lwd_axis)
for (i in seq_along(arms))
  jrnl_series(h0[[i]]$T, h0[[i]]$mean, arms[[i]]$col, lty = arms[[i]]$lty)
jrnl_axes(at_x = xt, at_y = yt_a)
jrnl_label_end(h0[[1]]$T, h0[[1]]$mean, "local", JCOL$local, pos = 3, off = 0.15)
jrnl_label_end(h0[[3]]$T, h0[[3]]$mean, "pMOM",  JCOL$mom,   pos = 3, off = 0.15)
jrnl_label_end(h0[[4]]$T, h0[[4]]$mean, "piMOM", JCOL$imom,  pos = 3, off = 0.15)
jrnl_panel("a")

# ---- (b) evidence FOR a true break ------------------------------------------
# All four coincide and grow linearly in T: the non-local prior gives up nothing
# where a break is real. Without this panel, (a) invites "it is merely conservative".
h1 <- lapply(arms, function(a) get(a$key, "H1"))
ylim_b <- c(0, max(unlist(lapply(h1, function(d) d$mean))) * 1.08)

plot(NA, xlim = range(TS), ylim = ylim_b, xlab = "Time periods (T)",
     ylab = "log Bayes factor", main = "", axes = FALSE, log = "x")
jrnl_grid(h = seq(0, 120, 20))
for (i in seq_along(arms))
  jrnl_series(h1[[i]]$T, h1[[i]]$mean, arms[[i]]$col, lty = arms[[i]]$lty)
jrnl_axes(at_x = xt, at_y = seq(0, 120, 20))
jrnl_legend("topleft",
            legend = sapply(arms, `[[`, "lab"),
            col    = sapply(arms, `[[`, "col"),
            lty    = sapply(arms, `[[`, "lty"),
            pch    = 16)
jrnl_panel("b")

# ---- (c) evidence vs the multiplicity burden --------------------------------
# Both quantities are log-odds, so they share one axis (never a dual axis).
# K = Ni*(T-3) candidates, so multiplicity contributes ~ +log K of spurious
# pressure. A prior controls false positives only if it supplies more evidence
# than the burden. Local: 0.5*log T against log(10T) -- the log terms cancel and
# the burden wins at every T. piMOM: ~1.55*sqrt(T) -- dominates and diverges.
logK <- log(Ni * (TS - 3))
ev   <- lapply(h0, function(d) -d$mean)          # evidence supplied = |log BF|
ylim_c <- c(0, max(unlist(ev)) * 1.08)

plot(NA, xlim = range(TS), ylim = ylim_c, xlab = "Time periods (T)",
     ylab = "Evidence vs. burden (log-odds)", main = "", axes = FALSE, log = "x")
jrnl_grid(h = seq(0, 50, 10))
for (i in seq_along(arms))
  jrnl_series(TS, ev[[i]], arms[[i]]$col, lty = arms[[i]]$lty)
lines(TS, logK, col = JCOL$ink, lty = 2, lwd = JRNL$lwd_ref)
jrnl_axes(at_x = xt, at_y = seq(0, 50, 10))
# At the right edge the series are well separated (local 3.0, normal slab 6.9,
# log K 8.8, pMOM 11.2, piMOM 38.8), so labels go below/above to miss each other.
jrnl_label_end(TS, logK,    "log K", JCOL$ink,   pos = 1, off = 0.25)
jrnl_label_end(TS, ev[[3]], "pMOM",  JCOL$mom,   pos = 3, off = 0.15)
jrnl_label_end(TS, ev[[4]], "piMOM", JCOL$imom,  pos = 1, off = 0.25)
jrnl_label_end(TS, ev[[1]], "local", JCOL$local, pos = 1, off = 0.25)
jrnl_panel("c")

dev.off()
cat("wrote", OUT, "\n")

# ---- numbers for the text ---------------------------------------------------
cat("\nRate table (for the manuscript):\n")
print(format(z$rates, digits = 3), row.names = FALSE)
cat(sprintf("\ntau: iMOM %.4f | pMOM %.4f | normal slab %.1f\n",
            z$tau["imom"], z$tau["mom"], z$tau["normal"]))
bf <- sapply(arms, function(a) exp(tail(get(a$key, "H0")$mean, 1)))
cat(sprintf("BF against one spurious break at T = %d:\n", max(TS)))
for (i in seq_along(arms)) cat(sprintf("  %-20s %.3e\n", arms[[i]]$lab, bf[i]))
cat(sprintf("\npiMOM vs Zellner at T = %d: %.3g x less evidence for a spurious break\n",
            max(TS), bf[1] / bf[4]))
