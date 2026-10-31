# Installs the bespoke, stripped-down mombf used by estimate_bisam() (based on mombf 3.5.4).
# Note: this replaces any CRAN mombf installed in the same library. Scripts that need the
# full package (e.g. scripts/16_prior_comparison_rates.R) should be run with CRAN mombf.
try(detach("package:mombf", unload = TRUE), silent = TRUE)
try(remove.packages("mombf"), silent = TRUE)
# .rs.restartR() # can help

Rcpp::compileAttributes(pkgdir = "./inst/mombf_3.5.4/mombf")
remotes::install_local("./inst/mombf_3.5.4/mombf", force = TRUE)
library(mombf)
