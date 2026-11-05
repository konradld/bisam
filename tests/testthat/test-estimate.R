sim_panel <- function(seed, ...) {
  set.seed(seed)
  contr_sim_breaks(n = 4, t = 20, nx = 1, iis = FALSE, sis = TRUE, const = TRUE,
                   pos.outl = 0, pos.step = c(10, 52), outl.mean = 0, step.mean = 4,
                   error.sd = 0.5, ...)
}
fit_quiet <- function(...) { out <- NULL; utils::capture.output(out <- estimate_bisam(...)); out }

test_that("estimate_bisam runs, is reproducible with set.seed and finds the breaks", {
  sim <- sim_panel(1)
  args <- list(data = sim$data, do_constant = TRUE, Ndraw = 600L, Nburn = 100L, beta_prior = "f",
               step_size_prior = "imom", step_incl_prior = "bern",
               step_size_scale = priorp2g(0.01, 1, prior = "iMom"))
  set.seed(2); f1 <- do.call(fit_quiet, args)
  set.seed(2); f2 <- do.call(fit_quiet, args)
  expect_s3_class(f1, "ism")
  expect_identical(f1$draws$omega, f2$draws$omega)
  expect_identical(f1$draws$beta, f2$draws$beta)
  true_cols <- rownames(sim$tr.idx)
  expect_true(all(f1$coefs$omega[true_cols] > 0.5))
  expect_true(mean(f1$coefs$omega[setdiff(names(f1$coefs$omega), true_cols)]) < 0.1)
  expect_equal(nrow(f1$draws$beta), 500L)
})

test_that("sparse and dense computation give the same posterior", {
  sim <- sim_panel(3)
  run <- function(sp) { set.seed(5); fit_quiet(data = sim$data, do_constant = TRUE, Ndraw = 3000L, Nburn = 500L,
    beta_prior = "g", step_size_prior = "imom", step_incl_prior = "bern", step_size_scale = 3,
    do_sparse_computation = sp) }
  d <- run(FALSE); s <- run(TRUE)
  expect_equal(unname(d$coefs$beta), unname(s$coefs$beta), tolerance = 0.05)
  expect_equal(d$coefs$sigma2, s$coefs$sigma2, tolerance = 0.05)
  expect_lt(max(abs(d$coefs$omega - s$coefs$omega)), 0.1)
})

test_that("sparse least squares equals lm.fit", {
  set.seed(1); X <- cbind(1, matrix(rnorm(200), 50)); y <- rnorm(50)
  expect_equal(sparse_lm_coef(X, matrix(y)), unname(stats::lm.fit(X, y)$coefficients), tolerance = 1e-10)
})

test_that("inv_gamma_params matches the invgamma definition", {
  # invgamma::pinvgamma(q, shape, scale = b) = pgamma(1/q, shape, rate = 1/b, lower.tail = FALSE)
  for (x in c(0.3, 1, 4)) {
    r <- inv_gamma_params(shape = 3, x = x, p = 0.9)
    expect_equal(stats::pgamma(1 / x, 3, rate = 1 / r$scale, lower.tail = FALSE), 0.9, tolerance = 1e-6)
    expect_equal(r$rate, 1 / r$scale)
  }
})

test_that("pip_window works with its default op", {
  sim <- sim_panel(4)
  set.seed(1)
  f <- fit_quiet(data = sim$data, do_constant = TRUE, Ndraw = 400L, Nburn = 100L, beta_prior = "f",
                 step_size_prior = "imom", step_incl_prior = "bern", step_size_scale = 3)
  w <- utils::capture.output(res <- pip_window(f, win_size = 2, pip_threshold = 0.5))
  expect_s3_class(res, "data.frame")
  expect_true(all(res$PIP >= 0.5))
})

test_that("contr_sim_breaks keeps fixed effects without covariates", {
  set.seed(1)
  s <- suppressWarnings(contr_sim_breaks(n = 3, t = 10, nx = 0, iis = FALSE, sis = TRUE, const = FALSE,
                                         ife = TRUE, tfe = TRUE, pos.outl = 0, pos.step = 5,
                                         outl.mean = 0, step.mean = 3, error.sd = 0.01))
  y <- matrix(s$data[, "y"], 10)
  # unit fixed effects show up as differences in unit means
  expect_gt(stats::sd(colMeans(y)), 1)
  expect_length(s$true.ife, 2)  # first unit effect dropped for identification
})

test_that("pip_window handles units shorter than the window", {
  set.seed(4)
  nm <- c(paste0("sis.A.", 3:8), paste0("sis.B.", 3:5))   # unit B has 3 periods < win_size 4
  om <- matrix(rbinom(200 * length(nm), 1, 0.3), 200, dimnames = list(NULL, nm))
  mod <- list(draws = list(omega = om, sis = om * rnorm(length(om))))
  utils::capture.output(res <- pip_window(mod, win_size = 4, pip_threshold = 0))
  ws <- mod$draws$sis[, grep("^sis\\.B", nm)]
  s <- sign(mean(rowSums(ws)))
  expect_equal(res$purity[res$unit == "B"], sum(sign(colMeans(ws)) == s) / 3 * s)
})
