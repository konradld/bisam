# Implementation notes: the compiled core of bisam

The compiled routines in `src/bisam_core.cpp` are derived from
[mombf](https://github.com/davidrusi/mombf) 3.5.4 (David Rossell et al., GPL >= 2). Only the
parts needed by `estimate_bisam()` were kept; they were re-implemented for speed and extended.
Earlier versions of this repository shipped them as a separate package (first `mombf`
3.5.4.9000, then `bisamcore`).

Internal entry points (not exported; call as `bisam:::name`):

| Function | Purpose |
|---|---|
| `bisam_ms_gibbs(y, x, prior, tau, phi, priorDelta, p_incl, bb_a, bb_b, niter, burnin, deltaini, imom_method = 0)` | Gibbs sampler over the break indicators of one unit (same algorithm as `mombf::modelSelection(..., method = "ALA", phi = phi)`). |
| `bisam_rnlp(y, x, prior, tau, phi, niter, burnin, thinit, newflag)` | Break-size draws from the MOM / iMOM posterior with known `phi` (Rossell & Telesca 2017, Alg. 2). |
| `bisam_log_joint(...)`, `bisam_imom_mode(...)`, `bisam_imom_stats()`, `bisam_set_imom_exhaustive_max()` | Diagnostics. |

Codes: `prior` 0 = product MOM (r = 1), 1 = product iMOM; `priorDelta` 1 = Bernoulli, 2 =
Beta-Bernoulli; `imom_method` 0 = ALA at the legacy mombf mode (default), 1 = ALA at the
near-global mode, 2 = Laplace at the near-global mode.

## iMOM integrated likelihood and the posterior mode

For the iMOM prior the integrated likelihood is approximated at the posterior mode. Under the
collinear step-dummy design the iMOM posterior is strongly multimodal (each sign pattern of the
break sizes can hold its own local mode). mombf's iteration, started at a ridge estimate, often
stops at a poor local mode: in tests on BISAM-like data the objective at its mode was a median of
2-3 (tau = 3.32) to 30-40 (tau = 17.3) log-units worse than the best mode, up to several hundred
for larger models. The near-global search (`imom_method = 1, 2`) is deterministic and never worse
than the legacy mode: models with up to 6 breaks are solved exactly over all sign patterns, larger
models by a sign-flip search (matched an exhaustive multi-start reference in ~98-100% of test
models). Against the exact integrated likelihood (importance sampling over all modes, m <= 5),
the error of the log Bayes factor for adding one break was:

| approximation | tau = 3.32: mean [sd] | tau = 17.3: mean [sd] |
|---|---|---|
| ALA at legacy mode (`imom_method = 0`, default) | -4.3 [6.7] | -18.2 [21.6] |
| ALA at near-global mode (`1`) | -2.0 [0.44] | -2.1 [0.57] |
| Laplace at near-global mode (`2`) | -0.11 [0.33] | -0.07 [0.28] |

ALA's remaining bias (about -2 log-units per break, from the BIC-type term -0.5 log n per break)
acts as an extra sparsity penalty. With Laplace the approximation is close to exact, which can
give much denser posteriors under the same model-space prior, so priors calibrated with ALA may
need recalibration before switching.

## Differences to mombf 3.5.4

* Only Gaussian errors with known residual variance, MOM/iMOM priors and the ALA
  approximation are supported (everything BISAM uses). All other model families, priors,
  methods, GGMs, survival, GLMs, etc. were removed.
* iMOM posterior mode: `imom_method = 0` runs the same iteration as mombf (ridge start,
  Gauss-Seidel coordinate sweeps, smallest stationary point on the current sign side, tolerance
  1e-5, max 50 sweeps), with the quartic root found by a bracketed Newton search warm-started at
  the coordinate's current value instead of building a `Polynomial` and running Jenkins-Traub at
  every update. Log integrated likelihoods agree with mombf to < 2e-5 (the solver tolerance).
  `imom_method = 1, 2` use the near-global mode instead (see above); they are 1.5-3.5x slower.
* Random numbers come from R's RNG (`set.seed()` works). mombf 3.5.4 drew the model-selection
  and break-size random numbers from an internal generator with fixed seeds.
* Break-size sampler (`bisam_rnlp`): truncated-normal draws use log-scale tail probabilities and
  quantiles (mombf computed 1 - Phi in the upper tail and clamped quantiles at 8.2 sd: for an iMOM
  prior with tau = 100 and a long step this biased the posterior mean by 0.6 sd and the sd by
  -27%); the inverse iMOM penalty is solved to machine precision (mombf: 1e-5); breaks that just
  entered the model can start at their conditional posterior mode. The stationary distribution
  was checked against exact 1- and 2-d posteriors (numerical integration), and the tail sampler
  against exact truncated normals and Robert's (1995) exponential-rejection sampler up to 40 sd.
  The sampler mixes slowly for nearby (correlated) breaks. Because the break indicators are drawn
  with the sizes integrated out, `estimate_bisam()` now runs `rnlp_iter = 25` iterations when a
  unit's set of breaks changed and 1 iteration otherwise (previously always 3, starting new breaks
  from stale values, which in full runs with many breaks inflated sigma^2 by ~5%).

The algorithms were validated against mombf 3.5.4 (identical Gibbs draws given the same
random numbers); see the commit history of the BISAM repository.
