try(detach("package:mombf", unload=TRUE))
try(detach("package:cli", unload=TRUE))
try(remove.packages("mombf"))
# .rs.restartR() # can help

Rcpp::compileAttributes(pkgdir = "./inst/mombf_3.5.4/mombf")
remotes::install_local("./inst/mombf_3.5.4/mombf", force = TRUE)
library(mombf)
