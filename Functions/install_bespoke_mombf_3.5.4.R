try(detach("package:mombf", unload=TRUE))
try(detach("package:BISAM", unload=TRUE))
try(detach("package:cli", unload=TRUE))
try(remove.packages("mombf"))
# .rs.restartR()

Rcpp::compileAttributes(pkgdir = "./Functions/mombf_3.5.4/mombf")
remotes::install_local("./Functions/mombf_3.5.4/mombf", force = TRUE)
library(mombf)
