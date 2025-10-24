try(detach("package:mombf", unload=TRUE))
try(detach("package:BISAM", unload=TRUE))
try(detach("package:cli", unload=TRUE))
# remove.packages("cli")
try(remove.packages("mombf"))
# rm(list = ls())
.rs.restartR()

Rcpp::compileAttributes(pkgdir = "~/BBD/00_code/mombf_3.5.4/mombf")
remotes::install_local("~/BBD/00_code/mombf_3.5.4/mombf", force = TRUE)

library(mombf)







# Only needed if you added/removed [[Rcpp::export]] functions
Rcpp::compileAttributes("~/BBD/00_code/mombf_3.5.4/mombf")

# Fast reload - compiles only changed files and loads
# devtools::clean_dll("~/BBD/00_code/mombf_3.5.4/mombf")  # Remove compiled code
devtools::load_all("~/BBD/00_code/mombf_3.5.4/mombf")
