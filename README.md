# bisam: Bayesian Indicator Saturated Models

R package for Bayesian detection of structural breaks (step shifts) and outliers in panel
data. Break indicators are sampled with non-local (MOM / iMOM) priors on the break sizes;
posterior inclusion probabilities (PIPs) identify breaks, optionally aggregated over time
windows. The repository also contains the scripts for the simulation study (including a
comparison with GETS indicator saturation) and the emissions application.

## Installation

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("konradld/bisam", ref = "dev_speedup")
```

or from a local clone (run in the repository root):

```r
remotes::install_local(".", force = TRUE)
```

The package contains C++ code, so a compiler is needed (on Windows: Rtools). Restart R after
reinstalling. `scripts/check_bisam_speed.R` checks that the installed version matches the
repository and times `estimate_bisam()`.

Earlier versions of this repository installed the compiled core as a separate package
(`mombf` 3.5.4.9000, later `bisamcore`); these are no longer needed and can be removed with
`remove.packages(c("bisamcore"))`. A CRAN `mombf` is unaffected (and is still used by
`scripts/16_prior_comparison_rates.R`).

## Usage

```r
library(bisam)
set.seed(1)
sim <- contr_sim_breaks(n = 10, t = 30, nx = 0, iis = FALSE, sis = TRUE, const = TRUE,
                        pos.outl = 0, pos.step = c(43, 108), outl.mean = 0,
                        step.mean = 3, error.sd = 1)
fit <- estimate_bisam(sim$data, do_constant = TRUE, Ndraw = 5000, Nburn = 1000,
                      beta_prior = "f", step_size_prior = "imom", step_incl_prior = "bern",
                      step_size_scale = priorp2g(0.01, 1, prior = "iMom"))
breaks <- pip_window(fit, win_size = 1, op = ">=", pip_threshold = 0.5)
```

See `?estimate_bisam`, `?pip_window`, `?contr_sim_breaks` and `?priorp2g`.

## Repository structure

| Path | Content |
|---|---|
| `R/`, `src/`, `man/`, `tests/`, `DESCRIPTION`, `NAMESPACE` | The `bisam` package. |
| `inst/IMPLEMENTATION.md` | Notes on the compiled core: what it took from mombf 3.5.4, what changed, and how it was validated. |
| `scripts/` | Simulation study (`00_`-`17_`), emissions application (`21_`, `22_`) and cluster versions (`adapted_for_cluster/`). Not part of the package. |
| `data/` | Emissions data used by the application scripts. |
| `output/` | Results of the simulation study and the application. |

The compiled core is derived from [mombf](https://github.com/davidrusi/mombf) 3.5.4
(D. Rossell et al.), hence the GPL (>= 2) licence.

## References

tba

## Contact

<lucas.konrad@wu.ac.at>

## Citation

If you use this code in your research, please cite:

```
tba
```
