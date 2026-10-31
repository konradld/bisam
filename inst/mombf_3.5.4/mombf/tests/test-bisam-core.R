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
cat("All BISAM core tests passed\n")
