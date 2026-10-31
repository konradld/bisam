# mombf (bespoke BISAM build, based on mombf 3.5.4)

A stripped-down version of [mombf](https://github.com/davidrusi/mombf) 3.5.4 (David Rossell
et al., GPL >= 2) that contains only what `R/estimate_bisam_fun.R` needs. The package keeps
the name `mombf` so that the BISAM scripts (`library(mombf)`, `priorp2g()`) work unchanged.
**Installing it replaces any CRAN mombf in the same library.** Scripts that need the full
package (e.g. `scripts/16_prior_comparison_rates.R`, which calls `nlpMarginal()` with Zellner
and normal priors) should be run with CRAN mombf.

## Contents

| Function | Purpose |
|---|---|
| `bisam_ms_gibbs(y, x, prior, tau, phi, priorDelta, p_incl, bb_a, bb_b, niter, burnin, deltaini)` | Gibbs sampler over the inclusion indicators of one unit's break columns. Same algorithm and output as `mombf::modelSelection(..., method = "ALA", phi = phi, enumerate = FALSE, initSearch = "none", center = FALSE, scale = FALSE)` for `momprior()` / `imomprior()` and `modelbinomprior()` / `modelbbprior()`. Returns `postSample` (logical matrix) and `margpp`. |
| `bisam_rnlp(y, x, prior, tau, phi, niter, burnin, thinit)` | Break-size draws from the MOM / iMOM posterior with known `phi`; returns the last of `niter` draws. Same algorithm as the earlier bespoke `rnlp_new(..., knownphi = TRUE, use_thinit = TRUE)`. |
| `bisam_log_joint(y, x, sel, prior, tau, phi, priorDelta, p_incl, bb_a, bb_b)` | Log integrated likelihood + log model prior of one model (diagnostics). |
| `priorp2g`, `mode2g`, `pmom`, `qmom`, `pimom`, `qimom` | Prior elicitation helpers, verbatim from mombf. |

Codes: `prior` 0 = product MOM (r = 1), 1 = product iMOM; `priorDelta` 1 = Bernoulli(`p_incl`),
2 = Beta-Bernoulli(`bb_a`, `bb_b`).

## Differences to mombf 3.5.4

* Only Gaussian errors with known residual variance, MOM/iMOM priors and the ALA
  approximation are supported (everything BISAM uses). All other model families, priors,
  methods, GGMs, survival, GLMs, etc. were removed.
* iMOM posterior mode: same iteration as mombf (ridge start, Gauss-Seidel coordinate sweeps,
  smallest stationary point on the current sign side, tolerance 1e-5, max 50 sweeps), but the
  quartic root is found by a bracketed Newton search warm-started at the coordinate's current
  value instead of building a `Polynomial` and running Jenkins-Traub at every update. Log
  integrated likelihoods agree with mombf to < 2e-5 (the solver tolerance).
* Random numbers come from R's RNG (`set.seed()` works). mombf 3.5.4 drew the model-selection
  and break-size random numbers from an internal generator with fixed seeds.
* Normal cdf/quantiles in the truncated-normal sampler use R's `pnorm`/`qnorm`.

The algorithms were validated against mombf 3.5.4 (identical Gibbs draws given the same
random numbers); see the commit history of the BISAM repository.
