# Basic checks of the BISAM core (run with: Rscript tests/test-bisam-core.R)
library(mombf)
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
cat("All BISAM core tests passed\n")
