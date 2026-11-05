# Basic checks of the BISAM core (run with: Rscript tests/test-bisam-core.R)
library(bisamcore)
set.seed(1)
t <- 30
z <- lower.tri(matrix(1, t, t), diag = TRUE)[, -c(1:2, t)] * 1.0
K <- ncol(z)
y <- as.vector(z[, c(5, 18)] %*% c(3, -2) + rnorm(t))

# MOM ALA log integrated likelihood against a direct R computation (phi = 1, tau = 2)
sel <- c(5, 9, 18); tau <- 2; x <- z[, sel]
S <- crossprod(x) + diag(1 / tau, 3); m <- solve(S, crossprod(x, y))
ref <- sum(log(m^2)) - 0.5 * (sum(y^2) - sum(m * (S %*% m))) -
  (0.5 * (t * log(2 * pi) + determinant(S)$modulus + 3 * log(tau)) + 3 * log(tau)) +
  3 * log(0.5) + (K - 3) * log(0.5)
stopifnot(abs(bisam_log_joint(y, z, sel, 0L, tau, 1, 1L, 0.5, 1, 1) - ref) < 1e-8)

# iMOM: the returned value is finite and the empty model is the Gaussian likelihood
tau_i <- priorp2g(0.05, 3, nu = 1, prior = "iMom")
stopifnot(is.finite(bisam_log_joint(y, z, sel, 1L, tau_i, 1, 1L, 0.5, 1, 1)))
stopifnot(abs(bisam_log_joint(y, z, integer(0), 1L, tau_i, 1, 2L, 0.5, 1, 1) -
                (sum(dnorm(y, 0, 1, log = TRUE)) + lbeta(1, 1 + K) - lbeta(1, 1))) < 1e-10)

# Gibbs output and reproducibility under set.seed()
run <- function() { set.seed(7); bisam_ms_gibbs(y, z, 1L, tau_i, 1, 1L, 0.5, 1, 1, 2L, 1L, rep(FALSE, K)) }
a <- run(); b <- run()
stopifnot(identical(a, b), is.logical(a$postSample), all(dim(a$postSample) == c(1, K)),
          all(a$margpp >= 0 & a$margpp <= 1), a$margpp[5] > 0.5, a$margpp[18] > 0.5)

# Break-size draws: correct length, non-zero, reproducible
r1 <- { set.seed(3); bisam_rnlp(y, z[, c(5, 18)], 1L, tau_i, 1, 3L, 2L, c(1, -1)) }
r2 <- { set.seed(3); bisam_rnlp(y, z[, c(5, 18)], 1L, tau_i, 1, 3L, 2L, c(1, -1)) }
stopifnot(length(r1) == 2, all(r1 != 0), identical(r1, r2))
# iMOM posterior mode: near-global search vs legacy and vs brute force over all orthants
set.seed(11)
tau3 <- priorp2g(0.01, 1, nu = 1, prior = "iMom")
for (rep in 1:25) {
  yy <- as.vector(z[, sample(K, 2)] %*% rnorm(2, 0, 3) + rnorm(t)); m <- sample(1:5, 1); s0 <- sort(sample(K, m))
  leg <- bisam_imom_mode(yy, z, s0, tau3, 1, 0L, numeric(0)); glo <- bisam_imom_mode(yy, z, s0, tau3, 1, 1L, numeric(0))
  stopifnot(glo$f <= leg$f + 1e-9)
  S <- as.matrix(expand.grid(rep(list(c(-1, 1)), m))); ref <- Inf
  for (r in seq_len(nrow(S))) for (mag in c(sqrt(tau3), 3 * sqrt(tau3)))
    ref <- min(ref, bisam_imom_mode(yy, z, s0, tau3, 1, 2L, S[r, ] * mag)$f)
  stopifnot(glo$f <= ref + 1e-8)
  A <- crossprod(z[, s0, drop = FALSE]); b <- as.vector(crossprod(z[, s0, drop = FALSE], yy)); th <- glo$th
  grad <- as.vector(A %*% th - b) - 2 * tau3 / th^3 + 2 / th
  H <- A; diag(H) <- diag(H) + 6 * tau3 / th^4 - 2 / th^2
  stopifnot(max(abs(grad)) < 1e-6, min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) > 0)
  # ALA and Laplace log joints from the mode
  k <- 0.5 * (m * log(tau3) - sum(yy^2) - t * log(2 * pi) - m * log(pi)) + K * log(0.5)
  stopifnot(abs(bisam_log_joint(yy, z, s0, 1L, tau3, 1, 1L, 0.5, 1, 1, 1L) - (k - glo$f - 0.5 * m * log(t))) < 1e-8,
            abs(bisam_log_joint(yy, z, s0, 1L, tau3, 1, 1L, 0.5, 1, 1, 0L) - (k - leg$f - 0.5 * m * log(t))) < 1e-8,
            abs(bisam_log_joint(yy, z, s0, 1L, tau3, 1, 1L, 0.5, 1, 1, 2L) -
                (k - glo$f + m / 2 * log(2 * pi) - 0.5 * as.numeric(determinant(H)$modulus))) < 1e-6)
}
# Break-size sampler (bisam_rnlp): stationary distribution vs the exact 1-d posterior,
# including a case far in the Normal tail (iMOM, tau = 100), where mombf's sampler was biased
lpr <- list(mom  = function(th, tau) log(th^2 / tau) + dnorm(th, 0, sqrt(tau), log = TRUE),
            imom = function(th, tau) 0.5 * log(tau) - 0.5 * log(pi) - log(th^2) - tau / th^2)
for (cs in list(list("mom", 0.348, 25, 0.6), list("imom", 3.317448, 25, 0), list("imom", 100, 28, 0.2))) {
  set.seed(4); pr <- cs[[1]]; tau1 <- cs[[2]]; k <- cs[[3]]
  x1 <- matrix(1, k, 1); y1 <- rnorm(k); y1 <- y1 - mean(y1) + cs[[4]]   # data mean fixed at cs[[4]]
  lf <- function(th) lpr[[pr]](th, tau1) - 0.5 * k * (th - mean(y1))^2
  g <- seq(-6, 6, length.out = 400001); g <- g[g != 0]; w <- exp(lf(g) - max(lf(g))); w <- w / sum(w)
  ex_m <- sum(w * g); ex_sd <- sqrt(sum(w * g^2) - ex_m^2)
  th <- mean(y1) + 0.5; dr <- numeric(20000)
  for (i in seq_along(dr)) { th <- bisam_rnlp(y1, x1, if (pr == "imom") 1L else 0L, tau1, 1, 1L, 0L, th); dr[i] <- th }
  dr <- dr[-(1:1000)]
  stopifnot(abs(mean(dr) - ex_m) < 0.1 * ex_sd, abs(sd(dr) / ex_sd - 1) < 0.1)
}

# Smart start (newflag): the flagged coordinate is set to the maximiser of its conditional posterior
x2 <- z[, c(5, 18)]; y2 <- as.vector(x2 %*% c(2, -1.5) + rnorm(t)); tau2 <- 3.317448
S2 <- crossprod(x2) + diag(1 / tau2, 2); xy2 <- as.vector(crossprod(x2, y2)); th1 <- 1.7
for (pr in c("imom", "mom")) {
  nlp_pr <- function(th2) { th <- c(th1, th2); 0.5 * sum(th * (S2 %*% th)) - sum(th * xy2) -
    sum(lpr[[pr]](th, tau2) - dnorm(th, 0, sqrt(tau2), log = TRUE)) }
  o <- list(optimize(nlp_pr, c(-10, -1e-4), tol = 1e-12), optimize(nlp_pr, c(1e-4, 10), tol = 1e-12))
  argbest <- o[[which.min(sapply(o, `[[`, "objective"))]]$minimum
  s0 <- bisam_rnlp(y2, x2, if (pr == "imom") 1L else 0L, tau2, 1, 0L, 0L, c(th1, 1), c(FALSE, TRUE))
  stopifnot(s0[1] == th1, abs(s0[2] - argbest) < 1e-5)
}
cat("All BISAM core tests passed\n")
