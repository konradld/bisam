# bisam 0.2.0

* The repository is now an R package (`library(bisam)`); the compiled core, previously the
  separate package `bisamcore` (and before that a bespoke `mombf` 3.5.4.9000), is part of it.
  Scripts load the package instead of sourcing `R/*.R`.
* No longer depends on `invgamma`, `nleqslv`, `MatrixModels`, `dplyr`, `extraDistr` or
  `RcppEigen`; `glmnet` (only for `beta_prior = "flasso"`) and `stochvol` (only for
  `sv_backend = "stochvol"`) are optional.
* Fixed: `pip_window()` with the default `op` failed in R >= 4.2; its handling of units with
  fewer periods than the window used an undefined index.
* Fixed: `contr_sim_breaks()` with `nx = 0`, both fixed effects and no constant dropped all
  fixed-effect columns, so the simulated data contained no fixed effects although
  `true.ife`/`true.tfe` were returned.
* Results of `estimate_bisam()` are bit-identical to the previous (sourced) version given the
  same seed.
