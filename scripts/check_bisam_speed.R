# ==============================================================================
# Quick check that the installed bisam package is current, and how long
# estimate_bisam() itself takes (without GETS or anything else).
# Run from the repository root:  Rscript scripts/check_bisam_speed.R
# Reference for exactly this setup (Linux, R 4.3, reference BLAS):
#   original mombf 3.5.4-based code: ~59 s;  current code: ~4 s.
# If you see times close to the first number, an old version is being used
# (reinstall the package and restart R).
# ==============================================================================
rm(list = ls())
if (!requireNamespace("bisam", quietly = TRUE)) {
  stop("bisam is not installed. From the repository root run: remotes::install_local('.')")
}
suppressMessages(library(bisam))
cat("bisam version:", as.character(packageVersion("bisam")), "\n")
cat("loaded from:  ", find.package("bisam"), "\n")
if (file.exists("DESCRIPTION") &&
    as.character(packageVersion("bisam")) != read.dcf("DESCRIPTION", "Version")[1]) {
  warning("The installed bisam differs from the version in this repository. Reinstall it and restart R.")
}

set.seed(1)
Ni <- 10; Nt <- 30; N_STEPS <- 1:4
POS_STEP_IN_Z <- sapply(N_STEPS, \(x) sample(1:(Nt - 3) + (x - 1) * (Nt - 3), 1))
POS_STEP <- POS_STEP_IN_Z + 2 * N_STEPS + (N_STEPS - 1)
sim <- contr_sim_breaks(n = Ni, t = Nt, nx = 0, iis = FALSE, sis = TRUE, const = TRUE,
                        pos.outl = 0, pos.step = POS_STEP, outl.mean = 0, step.mean = 3, error.sd = 1)

NDRAW <- 3000L
Rprof(prof_file <- tempfile(), interval = 0.01)
invisible(capture.output(res <- estimate_bisam(
  data = sim$data, do_constant = TRUE, Ndraw = NDRAW, Nburn = NDRAW %/% 5,
  beta_prior = "f", step_size_prior = "imom", step_incl_prior = "bern",
  beta_variance_scale = 100, step_incl_prob = 0.5,
  step_size_scale = priorp2g(0.01, 1, nu = 1, prior = "iMom"), do_split_Z = TRUE)))
Rprof(NULL)

secs <- as.numeric(res$meta$timer, units = "secs")
cat(sprintf("\nestimate_bisam(): %d draws in %.1f s (%.2f s per 1000 draws)\n", NDRAW, secs, 1000 * secs / NDRAW))
cat("true breaks:", POS_STEP_IN_Z, " | top PIPs:", head(order(-res$coefs$omega), 4), "\n\n")
cat("Where the time goes (top functions by self time):\n")
print(head(summaryRprof(prof_file)$by.self[, c("self.time", "self.pct")], 10))
