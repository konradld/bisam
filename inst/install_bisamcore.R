# Installs bisamcore, the compiled core used by estimate_bisam(). bisamcore is derived from
# mombf 3.5.4 (D. Rossell et al., GPL >= 2), stripped down to what BISAM needs and extended
# (see inst/bisamcore/README.md). It does not conflict with CRAN mombf, which can be
# installed alongside (needed e.g. by scripts/16_prior_comparison_rates.R).
try(detach("package:bisamcore", unload = TRUE), silent = TRUE)
try(remove.packages("bisamcore"), silent = TRUE)
# Earlier versions of this repository installed the same code under the name "mombf"
# (version 3.5.4.9000). Remove that build if present; a CRAN mombf is left untouched.
old_mombf <- tryCatch(as.character(utils::packageVersion("mombf")), error = function(e) "")
if (old_mombf == "3.5.4.9000") {
  try(detach("package:mombf", unload = TRUE), silent = TRUE)
  try(unloadNamespace("mombf"), silent = TRUE)
  try(remove.packages("mombf", lib = dirname(find.package("mombf"))), silent = TRUE)
  message("Removed the old bespoke 'mombf' (3.5.4.9000); it is replaced by 'bisamcore'.")
}
# .rs.restartR() # can help

Rcpp::compileAttributes(pkgdir = "./inst/bisamcore")
remotes::install_local("./inst/bisamcore", force = TRUE)
library(bisamcore)
