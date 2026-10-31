// Lean BISAM core: the only compiled code BISAM needs (bespoke, stripped-down mombf).
//
//  * bisam_ms_gibbs(): Gibbs sampler over the break-inclusion indicators of one unit,
//    Gaussian errors with KNOWN residual variance phi, product MOM (r = 1) or product
//    iMOM prior on the break sizes, Bernoulli or Beta-Bernoulli model prior. Same algorithm
//    as mombf 3.5.4 modelSelection(method = 'ALA', phi = ...). Integrated likelihoods:
//      - pMOM : ALA, plug-in of the Gaussian posterior mean into the MOM penalty
//               (mombf pmomMarginalKC, method 2)
//      - piMOM: posterior mode th* of the iMOM posterior, then (imom_method)
//          0  ALA at the mode reached by mombf's own iteration (legacy; often a poor
//             local mode under the collinear step design)
//          1  ALA at a near-global mode (default)
//          2  Laplace (Hessian) at the near-global mode
//        where ALA = log integrand at th* - m/2 log(n) (mombf pimomMarginalKC, method 2).
//    Legacy mode: mombf's iteration (ridge start, Gauss-Seidel sweeps, smallest same-sign
//    root of a quartic per coordinate) is reproduced exactly, but each root comes from a
//    bracketed Newton search warm-started at the coordinate's current value instead of a
//    Polynomial object + Jenkins-Traub; see imom_mode_legacy().
//    Near-global mode: see ImomSolver::solve(); deterministic, never worse than legacy.
//    Integrated likelihoods of visited models are cached within a call.
//
//  * bisam_rnlp(): break-size draws from the non-local posterior given the selected
//    model and known phi (port of mombf's rnlpPost_lm_new / rnlp_Gibbs).
//
//  * bisam_log_joint(), bisam_imom_mode(), bisam_imom_stats(),
//    bisam_set_imom_exhaustive_max(): diagnostics and tests.
//
// All random numbers come from R's RNG, so results are reproducible with set.seed().
// (mombf 3.5.4 used an internal generator with fixed seeds that set.seed() did not touch.)

#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
#include <chrono>

using namespace Rcpp;

static inline double bisam_unif() { return unif_rand(); }

static const double kLog2Pi = 1.8378770664093453;
static const double kLogPi  = 1.1447298858494002;

// ---------------------------------------------------------------------------
// Problem data for one unit
// ---------------------------------------------------------------------------
struct MsProblem {
  int n, K;            // observations, candidate break columns
  arma::mat XtX;       // K x K
  arma::vec ytX;       // K
  double sumy2, logly0; // sum(y^2) and log N(y; 0, phi I)
  double phi, tau;
  int prior;           // 0 = pMOM, 1 = piMOM
  int r;               // MOM power (fixed to 1 in BISAM)
  int priorDelta;      // 1 = binomial, 2 = beta-binomial
  double p_incl, bb_a, bb_b;
  int imom_method;     // piMOM: 0 = ALA at legacy mombf mode, 1 = ALA at near-global mode,
                       //        2 = Laplace (Hessian) at near-global mode
};

static double log_model_prior(const MsProblem &P, int k) {
  if (P.priorDelta == 1) {
    return k * std::log(P.p_incl) + (P.K - k) * std::log(1.0 - P.p_incl);
  }
  return R::lbeta(P.bb_a + k, P.bb_b + P.K - k) - R::lbeta(P.bb_a, P.bb_b);
}

// ---------------------------------------------------------------------------
// piMOM posterior mode (known phi)
//
// Objective (negative log posterior up to constants), th in R^m, all th_i != 0:
//   f(th) = .5 (th'A th - 2 b'th)/phi + c sum 1/th_i^2 + sum log th_i^2,  c = tau*phi
// ---------------------------------------------------------------------------
struct ImomObj {
  const arma::mat &A; const arma::vec &b; double phi, c;
};

static double imom_f(const ImomObj &O, const arma::vec &th) {
  double q = arma::dot(th, O.A * th), lin = arma::dot(O.b, th), pen = 0.0;
  for (arma::uword i = 0; i < th.n_elem; i++) {
    double t2 = th[i] * th[i];
    pen += O.c / t2 + std::log(t2);
  }
  return 0.5 * (q - 2.0 * lin) / O.phi + pen;
}

// Root search for one coordinate. With y = |th_i| on the side sign(th_i), stationary points of
// the 1-d section solve h(y) = a y^4 - B y^3 + 2 y^2 - 2c = 0 (B = sign * partial residual).
// h(0) = -2c < 0 and h'(y) = y (4a y^2 - 3B y + 4) has at most two positive roots yc1 < yc2,
// which bracket the roots of h. mombf took the first same-sign real root returned by a
// Jenkins-Traub solver, i.e. the smallest one; we return exactly that root, found by
// safeguarded Newton warm-started at the current |th_i| (no polynomial object, no allocation).
static double safeguarded_root(double a, double B, double c, double lo, double hi, double y) {
  if (!(y > lo && y < hi)) {  // warm start outside the bracket: start at the regula-falsi point
    double hlo = ((a * lo - B) * lo + 2.0) * lo * lo - 2.0 * c;
    double hhi = ((a * hi - B) * hi + 2.0) * hi * hi - 2.0 * c;
    y = (hhi > hlo) ? lo + (hi - lo) * (-hlo) / (hhi - hlo) : 0.5 * (lo + hi);
    if (!(y > lo && y < hi)) y = 0.5 * (lo + hi);
  }
  for (int it = 0; it < 200; it++) {
    double y2 = y * y;
    double h  = ((a * y - B) * y + 2.0) * y2 - 2.0 * c;
    if (h == 0.0) return y;
    if (h < 0.0) lo = y; else hi = y;
    double hp = y * ((4.0 * a * y - 3.0 * B) * y + 4.0);
    double ynew = (hp > 0.0) ? y - h / hp : 0.5 * (lo + hi);
    if (!(ynew > lo && ynew < hi)) ynew = 0.5 * (lo + hi);
    if (std::fabs(ynew - y) <= 4e-16 * y || hi - lo <= 4e-16 * hi) return ynew;
    y = ynew;
  }
  return y;
}

static double smallest_positive_root(double a, double B, double c, double ystart) {
  double hi = 1.0 + std::max(std::fabs(B), std::max(2.0, 2.0 * c)) / a; // Cauchy bound on |roots|
  double lo = 0.0;
  double disc = 9.0 * B * B - 64.0 * a;
  if (B > 0.0 && disc > 0.0) {
    double sq = std::sqrt(disc);
    double yc1 = (3.0 * B - sq) / (8.0 * a), yc2 = (3.0 * B + sq) / (8.0 * a);
    double h1 = ((a * yc1 - B) * yc1 + 2.0) * yc1 * yc1 - 2.0 * c;
    if (h1 > 0.0) hi = yc1; else lo = yc2;  // h < 0 on (0, yc2] when h(yc1) <= 0
  }
  return safeguarded_root(a, B, c, lo, hi, ystart);
}

// piMOM posterior mode, replicating mombf's imomModeK + imomIntegralApproxC exactly:
// ridge start (A + tau I)^{-1} b, Gauss-Seidel sweeps over coordinates, each coordinate set
// to the smallest stationary point on its current sign side, stop when the largest change in
// a sweep is <= 1e-5 or after 50 sweeps. Returns f at the final point.
// NOTE: the iMOM posterior under the collinear step design is strongly multimodal, so the
// starting point determines which local mode is found. Warm starting from another model's or
// iteration's mode would make the integrated likelihood depend on the chain's history; the
// ridge start keeps it a deterministic function of (data, model), as in mombf.
static double imom_mode_legacy(const ImomObj &O, double tau, arma::vec &th) {
  const arma::uword m = O.b.n_elem;
  arma::mat V = O.A;
  V.diag() += tau;
  arma::mat L;
  if (arma::chol(L, V, "lower")) {
    th = arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), O.b));
  } else {
    th = arma::solve(V, O.b);
  }
  arma::vec r = O.A * th;   // r = A th, kept current so each partial residual is O(1)
  double err = 1.0;
  for (int sweep = 0; sweep < 50 && err > 1e-5; sweep++) {
    err = 0.0;
    for (arma::uword i = 0; i < m; i++) {
      double ti = th[i];
      if (ti == 0.0) continue;  // mombf never moves an exact zero (no same-sign root)
      double a = O.A(i, i) / O.phi;
      double bb = (O.b[i] - (r[i] - O.A(i, i) * ti)) / O.phi;
      double s = (ti > 0.0) ? 1.0 : -1.0;
      double y = smallest_positive_root(a, s * bb, O.c, std::fabs(ti));
      double tnew = s * y;
      err = std::max(err, std::fabs(ti - tnew));
      r += O.A.col(i) * (tnew - ti);
      th[i] = tnew;
    }
  }
  return imom_f(O, th);
}

// ---------------------------------------------------------------------------
// Near-global piMOM posterior mode (imom_method 1 and 2)
//
// The legacy iteration stops at whatever local mode its ridge start leads to. Each sign
// pattern (orthant) of th can hold its own local mode, and under the strongly collinear step
// design the legacy one is frequently far from the best (tens to hundreds of log-units in
// tests). ImomSolver::solve() returns a mode that is never worse than the legacy one and is a
// deterministic function of (A, b, phi, tau), so the integrated likelihood remains a fixed
// function of (data, model) and the Gibbs sampler keeps a well-defined target:
//   1. legacy mode, polished to the exact local minimum (Newton, signs fixed);
//   2a. m <= kImomExhaustiveMax: every orthant, visited in Gray-code order, each optimised
//       from its neighbour's optimum (exact over sign patterns in all tests);
//   2b. larger m: greedy and OLS starts with sign-switching coordinate descent, then a
//       search over sign flips of contiguous blocks of coordinates until none improves.
// Validation (BISAM-like data, tau in {0.348, 3.32, 17.3}): 2a matched an exhaustive
// multi-start reference in all test models; 2b matched it in ~98-100%, worst gap < 1.
// ---------------------------------------------------------------------------
static inline double g_section(double a, double B, double c, double y) {
  return 0.5 * a * y * y - B * y + c / (y * y) + 2.0 * std::log(y);
}

// Global minimiser of g(y) = .5 a y^2 - B y + c/y^2 + 2 log y over y > 0; stores g in *gval.
// Local minima are the up-crossings of h(y) = a y^4 - B y^3 + 2 y^2 - 2c (h = y^3 g'):
// one root if h is increasing on y > 0, otherwise up to three roots split by the critical
// points yc1 < yc2 of h (up-crossings in (0, yc1) and (yc2, inf)).
static double best_positive_min(double a, double B, double c, double ystart, double *gval) {
  double hi = 1.0 + std::max(std::fabs(B), std::max(2.0, 2.0 * c)) / a; // Cauchy bound
  double disc = 9.0 * B * B - 64.0 * a;
  if (B <= 0.0 || disc <= 0.0) {
    double y = safeguarded_root(a, B, c, 0.0, hi, ystart);
    *gval = g_section(a, B, c, y);
    return y;
  }
  double sq = std::sqrt(disc);
  double yc1 = (3.0 * B - sq) / (8.0 * a), yc2 = (3.0 * B + sq) / (8.0 * a);
  double h1 = ((a * yc1 - B) * yc1 + 2.0) * yc1 * yc1 - 2.0 * c;
  double h2 = ((a * yc2 - B) * yc2 + 2.0) * yc2 * yc2 - 2.0 * c;
  double best_y = 0.0, best_g = R_PosInf;
  if (h1 > 0.0) {
    double y = safeguarded_root(a, B, c, 0.0, yc1, ystart), g = g_section(a, B, c, y);
    if (g < best_g) { best_g = g; best_y = y; }
  }
  if (h2 < 0.0 || h1 <= 0.0) {
    double y = safeguarded_root(a, B, c, yc2, std::max(hi, 2.0 * yc2), ystart), g = g_section(a, B, c, y);
    if (g < best_g) { best_g = g; best_y = y; }
  }
  *gval = best_g;
  return best_y;
}

// Small dense Cholesky (lower, in place, row-major m x m) and solve; false if not PD.
static bool chol_inplace(std::vector<double> &M, int m) {
  for (int j = 0; j < m; j++) {
    double d = M[j * m + j];
    for (int k = 0; k < j; k++) d -= M[j * m + k] * M[j * m + k];
    if (!(d > 0.0)) return false;
    d = std::sqrt(d);
    M[j * m + j] = d;
    for (int i = j + 1; i < m; i++) {
      double v = M[i * m + j];
      for (int k = 0; k < j; k++) v -= M[i * m + k] * M[j * m + k];
      M[i * m + j] = v / d;
    }
  }
  return true;
}
static void chol_solve(const std::vector<double> &L, int m, const double *rhs, double *x) {
  for (int i = 0; i < m; i++) {
    double v = rhs[i];
    for (int k = 0; k < i; k++) v -= L[i * m + k] * x[k];
    x[i] = v / L[i * m + i];
  }
  for (int i = m - 1; i >= 0; i--) {
    double v = x[i];
    for (int k = i + 1; k < m; k++) v -= L[k * m + i] * x[k];
    x[i] = v / L[i * m + i];
  }
}

// ---------------------------------------------------------------------------
// Allocation-free solver for the near-global mode. A and b are stored pre-divided by phi:
//   f(th) = .5 th'A th - b'th + sum_i (c/th_i^2 + log th_i^2)
// ---------------------------------------------------------------------------
static const int kCdSweeps = 30;       // free descent: coordinate sweeps before each Newton polish
static const int kCdSweepsFixed = 3;   // signs-fixed descent: Newton converges from there
static const double kCdTol = 1e-7;
static int kImomExhaustiveMax = 6;     // models up to this size: visit every orthant
static int kImomExStarts = 1;          // starts per orthant in the exhaustive search (1 or 2)

class ImomSolver {
public:
  int m; double c, sc;
  std::vector<double> A, H, L;                      // m x m
  std::vector<double> vecbuf;                       // all length-m work vectors, one allocation
  double *b, *r, *g, *step, *tn, *w1, *w2, *w3, *ols;

  ImomSolver(const arma::mat &A0, const arma::vec &b0, double phi, double c0)
    : m((int) b0.n_elem), c(c0), sc(std::sqrt(c0)), A(m * m), H(m * m), L(m * m), vecbuf(9 * m) {
    double *p = vecbuf.data();
    b = p; r = p + m; g = p + 2 * m; step = p + 3 * m; tn = p + 4 * m; w1 = p + 5 * m; w2 = p + 6 * m; w3 = p + 7 * m; ols = p + 8 * m;
    for (int i = 0; i < m; i++) { b[i] = b0[i] / phi; for (int j = 0; j < m; j++) A[i * m + j] = A0(i, j) / phi; }
    // (near-)OLS solution, for starting values (H used as scratch)
    double tr = 0.0; for (int i = 0; i < m; i++) tr += A[i * m + i];
    H = A;
    for (int i = 0; i < m; i++) H[i * m + i] += 1e-8 * (1.0 + tr / m);
    if (chol_inplace(H, m)) chol_solve(H, m, b, ols);
    else for (int i = 0; i < m; i++) ols[i] = b[i] / A[i * m + i];
    for (int i = 0; i < m; i++) if (ols[i] == 0.0) ols[i] = 1e-8;
  }

  double fval(const double *th) const {
    double q = 0.0, lin = 0.0, pen = 0.0;
    for (int i = 0; i < m; i++) {
      const double *Ai = &A[i * m];
      double s = 0.0;
      for (int j = 0; j < m; j++) s += Ai[j] * th[j];
      double t2 = th[i] * th[i];
      q += th[i] * s; lin += b[i] * th[i]; pen += c / t2 + std::log(t2);
    }
    return 0.5 * q - lin + pen;
  }
  void Ax(const double *th, double *out) const {
    for (int i = 0; i < m; i++) { const double *Ai = &A[i * m]; double s = 0.0; for (int j = 0; j < m; j++) s += Ai[j] * th[j]; out[i] = s; }
  }
  void add_col(double *out, int i, double d) const { for (int k = 0; k < m; k++) out[k] += A[k * m + i] * d; }

  // coordinate i: global (both signs) or own-sign update of the 1-d section; r = A th kept current
  double coord_update(double *th, int i, bool free_sign) {
    double a = A[i * m + i], bb = b[i] - (r[i] - a * th[i]);
    double ystart = std::max(std::fabs(th[i]), 1e-8), tnew, gp, gm;
    if (free_sign) {
      double yp = best_positive_min(a, bb, c, ystart, &gp), ym = best_positive_min(a, -bb, c, ystart, &gm);
      tnew = (gp <= gm) ? yp : -ym;
    } else {
      double s = (th[i] > 0.0) ? 1.0 : -1.0;
      tnew = s * best_positive_min(a, s * bb, c, ystart, &gp);
    }
    double d = tnew - th[i];
    add_col(r, i, d);
    th[i] = tnew;
    return std::fabs(d) / (1.0 + std::fabs(tnew));
  }

  // Newton / Levenberg-Marquardt with all signs fixed, to the exact local minimum.
  // f must equal fval(th); returns the final f (never larger).
  double newton(double *th, double f) {
    for (int it = 0; it < 200; it++) {
      Ax(th, r);
      double hmax = 0.0;
      for (int i = 0; i < m; i++) {
        double t = th[i], t2 = t * t;
        g[i] = r[i] - b[i] - 2.0 * c / (t2 * t) + 2.0 / t;
        for (int j = 0; j < m; j++) H[i * m + j] = A[i * m + j];
        H[i * m + i] += 6.0 * c / (t2 * t2) - 2.0 / t2;
        hmax = std::max(hmax, std::fabs(H[i * m + i]));
      }
      double lam = 0.0, fnew = f;
      bool accepted = false;
      for (int tries = 0; tries < 40 && !accepted; tries++) {
        L = H;
        if (lam > 0.0) for (int i = 0; i < m; i++) L[i * m + i] += lam;
        if (chol_inplace(L, m)) {
          chol_solve(L, m, g, step);
          double slope = 0.0;
          for (int i = 0; i < m; i++) slope += g[i] * step[i];
          double t = 1.0;
          for (int ls = 0; ls < 60; ls++) {
            bool samesign = true;
            for (int i = 0; i < m; i++) {
              tn[i] = th[i] - t * step[i];
              if (!(tn[i] * th[i] > 0.0)) { samesign = false; break; }
            }
            if (samesign) {
              fnew = fval(tn);
              if (fnew <= f - 1e-4 * t * slope) { accepted = true; break; }
            }
            t *= 0.5;
          }
        }
        if (!accepted) lam = (lam == 0.0) ? 1e-8 * (1.0 + hmax) : lam * 10.0;
      }
      if (!accepted) break;   // no descent possible at machine precision: local minimum
      double rel = 0.0;
      for (int i = 0; i < m; i++) { rel = std::max(rel, std::fabs(tn[i] - th[i]) / (1.0 + std::fabs(th[i]))); th[i] = tn[i]; }
      f = fnew;
      if (rel < 1e-13 && lam == 0.0) break;
    }
    return f;
  }

  // signs fixed: a few own-sign coordinate sweeps, then Newton
  double descend_fixed(double *th) {
    Ax(th, r);
    for (int sweep = 0; sweep < kCdSweepsFixed; sweep++) {
      double err = 0.0;
      for (int i = 0; i < m; i++) err = std::max(err, coord_update(th, i, false));
      if (err < kCdTol) break;
    }
    return newton(th, fval(th));
  }

  // free signs: global coordinate sweeps + Newton, repeated until a sweep at the polished
  // point no longer improves f (local minimum that is coordinatewise optimal over both signs)
  double descend_free(double *th) {
    double f = R_PosInf;
    for (int round = 0; round < 50; round++) {
      Ax(th, r);
      for (int sweep = 0; sweep < kCdSweeps; sweep++) {
        double err = 0.0;
        for (int i = 0; i < m; i++) err = std::max(err, coord_update(th, i, true));
        if (err < kCdTol) break;
      }
      f = newton(th, fval(th));
      std::copy(th, th + m, w3);
      Ax(w3, r);
      double e2 = 0.0;
      for (int i = 0; i < m; i++) e2 = std::max(e2, coord_update(w3, i, true));
      if (e2 < kCdTol || fval(w3) >= f) break;
      std::copy(w3, w3 + m, th);
    }
    return f;
  }

  // Near-global mode; `best` holds the legacy mode on input and the result on output.
  double solve(double *best) {
    double fbest = newton(best, fval(best));        // 1. legacy mode, polished
    auto better = [&](double f) { return f < fbest - 1e-10 * (1.0 + std::fabs(fbest)); };
    auto take = [&](const double *t, double f) { if (better(f)) { fbest = f; std::copy(t, t + m, best); return true; } return false; };
    if (m <= kImomExhaustiveMax) {
      // 3a. every orthant, in Gray-code order: consecutive orthants differ in one sign, and
      //     the descent in each starts from the previous orthant's optimum with that one
      //     coordinate flipped (warm start); optionally a second start at magnitudes
      //     |ols| + sqrt(c). The visiting order is fixed, so the result is deterministic.
      std::copy(best, best + m, w1);
      const unsigned long north = 1UL << m;
      for (unsigned long k = 1; k < north; k++) {
        int gi = __builtin_ctzl(k);
        w1[gi] = -w1[gi];
        double f1 = descend_fixed(w1);
        take(w1, f1);
        if (kImomExStarts > 1) {
          for (int i = 0; i < m; i++) w2[i] = (w1[i] > 0.0 ? 1.0 : -1.0) * (std::fabs(ols[i]) + sc);
          take(w2, descend_fixed(w2));
        }
      }
    } else {
      // 2. deterministic free starts: greedy sequential build-up from 0, and OLS
      std::fill(w1, w1 + m, 0.0); std::fill(r, r + m, 0.0);
      for (int i = 0; i < m; i++) coord_update(w1, i, true);
      take(w1, descend_free(w1));
      std::copy(ols, ols + m, w1);
      take(w1, descend_free(w1));
      // 3b. block sign-flip search: flip any contiguous block, best point of the new orthant
      //     with those signs locked (from the flipped point and from the flipped point with
      //     the block's magnitudes reset to sqrt(c)); accept strict improvements and release
      //     them to the free descent; repeat until no block flip improves.
      for (int pass = 0; pass < 4 * m + 10; pass++) {
        bool improved = false;
        for (int width = 1; width <= m; width++) {
          for (int i = 0; i + width <= m; i++) {
            std::copy(best, best + m, w1);
            for (int q = 0; q < width; q++) w1[i + q] = -w1[i + q];
            std::copy(w1, w1 + m, w2);
            for (int q = 0; q < width; q++) w2[i + q] = (w1[i + q] > 0.0 ? sc : -sc);
            double f1 = descend_fixed(w1), f2 = descend_fixed(w2);
            bool imp = (f1 <= f2) ? take(w1, f1) : take(w2, f2);
            if (imp) {
              improved = true;
              std::copy(best, best + m, w1);
              take(w1, descend_free(w1));
            }
          }
        }
        if (!improved) break;
      }
    }
    return fbest;
  }
};

static double imom_mode_global(const ImomObj &O, double tau, arma::vec &th) {
  th.set_size(O.b.n_elem);
  imom_mode_legacy(O, tau, th);
  ImomSolver S(O.A, O.b, O.phi, O.c);
  return S.solve(th.memptr());
}

// [[Rcpp::export]]
int bisam_set_imom_exhaustive_max(int mmax, int starts = -1) {
  // diagnostics / tuning: largest model solved exhaustively, and starts per orthant (1 or 2)
  int old = kImomExhaustiveMax;
  if (mmax < 0 || mmax > 20) stop("mmax must be in 0..20");
  kImomExhaustiveMax = mmax;
  if (starts == 1 || starts == 2) kImomExStarts = starts;
  return old;
}

// Diagnostics: piMOM mode of model `sel` (1-based) by method
//   0 = legacy (mombf), 1 = global, 2 = Newton from `start` with signs fixed,
//   3 = global coordinate descent + Newton from `start`
// [[Rcpp::export]]
List bisam_imom_mode(NumericVector y, NumericMatrix x, IntegerVector sel, double tau, double phi,
                     int method, NumericVector start) {
  int n = y.size(), K = x.ncol(), m = sel.size();
  arma::mat X(x.begin(), n, K, false, true);
  arma::vec Y(y.begin(), n, false, true);
  arma::uvec idx(m);
  for (int i = 0; i < m; i++) {
    if (sel[i] < 1 || sel[i] > K) stop("bisam_imom_mode: sel must index columns of x (1-based)");
    idx[i] = sel[i] - 1;
  }
  arma::mat Xs = X.cols(idx);
  arma::mat A = Xs.t() * Xs;
  arma::vec b = Xs.t() * Y;
  ImomObj O{A, b, phi, tau * phi};
  arma::vec th(m);
  double f;
  if (method == 0) {
    f = imom_mode_legacy(O, tau, th);
  } else if (method == 1) {
    f = imom_mode_global(O, tau, th);
  } else {
    if (start.size() != m) stop("bisam_imom_mode: start must have length(sel)");
    th = arma::vec(start.begin(), m);
    for (int i = 0; i < m; i++) if (th[i] == 0.0) stop("bisam_imom_mode: start must be non-zero");
    ImomSolver S(A, b, phi, tau * phi);
    f = (method == 2) ? S.newton(th.memptr(), S.fval(th.memptr())) : S.descend_free(th.memptr());
  }
  return List::create(_["th"] = NumericVector(th.begin(), th.end()), _["f"] = f);
}

// Diagnostics: number of iMOM mode searches and seconds spent, by model size
static double g_mode_calls[64], g_mode_secs[64];

// [[Rcpp::export]]
DataFrame bisam_imom_stats(bool reset = false) {
  IntegerVector m; NumericVector calls, secs;
  for (int i = 0; i < 64; i++) if (g_mode_calls[i] > 0) { m.push_back(i); calls.push_back(g_mode_calls[i]); secs.push_back(g_mode_secs[i]); }
  if (reset) for (int i = 0; i < 64; i++) g_mode_calls[i] = g_mode_secs[i] = 0.0;
  return DataFrame::create(_["m"] = m, _["calls"] = calls, _["seconds"] = secs);
}

// ---------------------------------------------------------------------------
// Log joint (log integrated likelihood + log model prior) for model `sel`
// ---------------------------------------------------------------------------
static double log_joint(const MsProblem &P, const std::vector<arma::uword> &sel) {
  const int m = (int) sel.size();
  double ans;
  if (m == 0) {
    ans = P.logly0;
  } else {
    arma::uvec idx(m);
    for (int i = 0; i < m; i++) idx[i] = sel[i];
    arma::mat A = P.XtX.submat(idx, idx);
    arma::vec b = P.ytX.elem(idx);
    if (P.prior == 1) { // piMOM, ALA = mode + BIC-type term (mombf: pimomMarginalKC, method 2)
      ImomObj O{A, b, P.phi, P.tau * P.phi};
      arma::vec th(m);
      auto t0 = std::chrono::steady_clock::now();
      double fopt = (P.imom_method == 0) ? imom_mode_legacy(O, P.tau, th) : imom_mode_global(O, P.tau, th);
      if (m < 64) {
        g_mode_calls[m] += 1.0;
        g_mode_secs[m] += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      }
      // log of the integrand at th is k - f(th):  N(y; X th, phi I) * prod_i iMOM(th_i; tau, phi)
      double k = 0.5 * (m * std::log(P.tau) - P.sumy2 / P.phi - P.n * kLog2Pi
                        - (P.n - m) * std::log(P.phi) - m * kLogPi);
      if (P.imom_method == 2) {
        // Laplace at the mode: log int exp(-f) ~ -f* + m/2 log(2 pi) - 1/2 log det H(th*)
        std::vector<double> H(m * m);
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) H[i * m + j] = A(i, j) / P.phi;
          double t2 = th[i] * th[i];
          H[i * m + i] += 6.0 * O.c / (t2 * t2) - 2.0 / t2;
        }
        if (!chol_inplace(H, m)) Rcpp::stop("bisam: iMOM Hessian at the mode is not positive definite");
        double logdetH = 0.0;
        for (int i = 0; i < m; i++) logdetH += 2.0 * std::log(H[i * m + i]);
        ans = k - fopt + 0.5 * m * kLog2Pi - 0.5 * logdetH;
      } else {
        // ALA as in mombf (method 2): mode + BIC-type term
        ans = k - fopt - 0.5 * m * std::log((double) P.n);
      }
    } else {            // pMOM, ALA = plug-in posterior mean (mombf: pmomMarginalKC, method 2)
      arma::mat S = A;
      S.diag() += 1.0 / P.tau;
      arma::mat L;
      if (!arma::chol(L, S, "lower")) Rcpp::stop("bisam_ms_gibbs: X'X + I/tau not positive definite");
      arma::vec mhat = arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), b));
      double logdetS = 2.0 * arma::accu(arma::log(L.diag()));
      double logtau = std::log(P.tau), logphi = std::log(P.phi), ldf = 0.0;
      for (int k = 3; k <= 2 * P.r - 1; k += 2) ldf += std::log((double) k);
      double num = -0.5 * (P.sumy2 - arma::dot(b, mhat)) / P.phi;
      double den = 0.5 * (P.n * (kLog2Pi + logphi) + logdetS + m * logtau)
                   + m * P.r * (logtau + logphi + ldf);
      double plug = 0.0;
      for (int i = 0; i < m; i++) plug += std::log(mhat[i] * mhat[i]);
      ans = P.r * plug + num - den;
    }
  }
  return ans + log_model_prior(P, m);
}

// ---------------------------------------------------------------------------
// Gibbs sampler over inclusion indicators (mirrors mombf::modelSelectionGibbs for
// family = normal, no groups / constraints / forced variables)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List bisam_ms_gibbs(NumericVector y, NumericMatrix x, int prior, double tau, double phi,
                    int priorDelta, double p_incl, double bb_a, double bb_b,
                    int niter, int burnin, LogicalVector deltaini, int imom_method = 1) {
  MsProblem P;
  P.n = y.size(); P.K = x.ncol();
  if (x.nrow() != P.n) stop("bisam_ms_gibbs: dimensions of y and x do not match");
  if (deltaini.size() != P.K) stop("bisam_ms_gibbs: deltaini must have length ncol(x)");
  arma::mat X(x.begin(), P.n, P.K, false, true);
  arma::vec Y(y.begin(), P.n, false, true);
  P.XtX = X.t() * X; P.ytX = X.t() * Y;
  P.sumy2 = arma::dot(Y, Y);
  P.logly0 = -0.5 * P.n * (kLog2Pi + std::log(phi)) - 0.5 * P.sumy2 / phi;
  P.phi = phi; P.tau = tau; P.prior = prior; P.r = 1;
  P.priorDelta = priorDelta; P.p_incl = p_incl; P.bb_a = bb_a; P.bb_b = bb_b;
  if (imom_method < 0 || imom_method > 2) stop("imom_method must be 0 (ALA, legacy mode), 1 (ALA, global mode) or 2 (Laplace, global mode)");
  P.imom_method = imom_method;

  const int K = P.K;
  std::string key(K, '0');
  std::vector<arma::uword> sel;
  for (int j = 0; j < K; j++) if (deltaini[j]) { sel.push_back(j); key[j] = '1'; }

  std::unordered_map<std::string, double> cache;  // log joint of every model visited
  double curJ = log_joint(P, sel);
  cache.emplace(key, curJ);

  const int ilow = -burnin, iupper = niter - burnin + 1;
  const int nsave = std::max(0, niter - burnin);
  LogicalMatrix postSample(nsave, K);
  NumericVector margpp(K);
  std::vector<arma::uword> selnew; selnew.reserve(K);
  int savecnt = 0;

  for (int i = ilow; i < iupper; i++) {
    for (int j = 0; j < K; j++) {
      bool in = key[j] == '1';
      selnew.clear();
      for (arma::uword s : sel) if ((int) s != j) selnew.push_back(s);
      if (!in) { // insert keeping ascending order
        selnew.insert(std::upper_bound(selnew.begin(), selnew.end(), (arma::uword) j), (arma::uword) j);
      }
      if ((int) selnew.size() > P.n) continue; // invalid model: no update, no draw
      key[j] = in ? '0' : '1';
      double newJ;
      auto it = cache.find(key);
      if (it != cache.end()) {
        newJ = it->second;
      } else {
        newJ = log_joint(P, selnew);
        cache.emplace(key, newJ);
      }
      double pp = std::exp(newJ - curJ);
      pp = pp / (1.0 + pp);
      if (i >= 0) margpp[j] += in ? (1.0 - pp) : pp;
      double u = bisam_unif();
      if (u < pp) {
        sel.swap(selnew); curJ = newJ;
      } else {
        key[j] = in ? '1' : '0';
      }
    }
    if (i > 0 && savecnt < nsave) {
      for (arma::uword s : sel) postSample(savecnt, s) = true;
      savecnt++;
    }
  }
  double denom = iupper - std::max(0, ilow);
  if (denom > 0) for (int j = 0; j < K; j++) margpp[j] /= denom;

  return List::create(_["postSample"] = postSample, _["margpp"] = margpp);
}

// Log integrated likelihood + log model prior of one model (testing / diagnostics)
// [[Rcpp::export]]
double bisam_log_joint(NumericVector y, NumericMatrix x, IntegerVector sel, int prior, double tau,
                       double phi, int priorDelta, double p_incl, double bb_a, double bb_b,
                       int imom_method = 1) {
  MsProblem P;
  P.n = y.size(); P.K = x.ncol();
  arma::mat X(x.begin(), P.n, P.K, false, true);
  arma::vec Y(y.begin(), P.n, false, true);
  P.XtX = X.t() * X; P.ytX = X.t() * Y; P.sumy2 = arma::dot(Y, Y);
  P.logly0 = -0.5 * P.n * (kLog2Pi + std::log(phi)) - 0.5 * P.sumy2 / phi;
  P.phi = phi; P.tau = tau; P.prior = prior; P.r = 1;
  P.priorDelta = priorDelta; P.p_incl = p_incl; P.bb_a = bb_a; P.bb_b = bb_b;
  if (imom_method < 0 || imom_method > 2) stop("imom_method must be 0 (ALA, legacy mode), 1 (ALA, global mode) or 2 (Laplace, global mode)");
  P.imom_method = imom_method;
  std::vector<arma::uword> selv;
  for (int s : sel) {
    if (s < 1 || s > P.K) stop("bisam_log_joint: sel must index columns of x (1-based)");
    selv.push_back(s - 1);
  }
  std::sort(selv.begin(), selv.end());
  return log_joint(P, selv);
}

// ---------------------------------------------------------------------------
// Break-size draws given the selected model (known phi)
//
// Port of mombf 3.5.4's rnlpPost_lm_new (knownphi branch) -> rnlp_Gibbs ->
// rtmvnormOutside_Gibbs -> rnorm_truncMult, for the product MOM (r = 1) and iMOM priors.
// Latent-variable Gibbs sampler: given th, draw a slice variable per coefficient from the
// prior penalty, which turns the non-local posterior into a N(m, phi S^{-1}) truncated to
// |th_j| >= upper_j; the truncated normal is then sampled one coordinate at a time in the
// whitened scale z = chol(S^{-1})^{-1} th / sqrt(phi). Matrix helpers keep mombf's loop order.
// The only change is that Normal cdf/quantiles come from R (pnorm/qnorm) instead of
// mombf's copy of dcdflib.
// ---------------------------------------------------------------------------
typedef std::vector<double> dvec;
typedef std::vector<dvec> dmat;   // row-major, dense p x p

static void cholesky_lower(const dmat &a, int n, dmat &L) {   // mombf choldc
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) L[i][j] = (j >= i) ? a[i][j] : 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      double sum = L[i][j];
      for (int k = i - 1; k >= 0; k--) sum -= L[i][k] * L[j][k];
      if (i == j) {
        L[i][i] = std::sqrt(sum);
      } else {
        L[j][i] = sum / std::max(std::fabs(L[i][i]), 1e-10);
      }
    }
  }
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) L[i][j] = 0.0;
}

static void invert_lower_inplace(dmat &L, int n) {              // mombf choldc_inv_internal
  for (int i = 0; i < n; i++) {
    L[i][i] = 1.0 / std::max(L[i][i], 1e-10);
    for (int j = i + 1; j < n; j++) {
      double sum = 0.0;
      for (int k = i; k < j; k++) sum -= L[j][k] * L[k][i];
      L[j][i] = sum / std::max(L[j][j], 1e-10);
    }
  }
}

static void inv_posdef(const dmat &a, int n, dmat &out) {      // mombf inv_posdef
  dmat b(n, dvec(n));
  cholesky_lower(a, n, b);
  invert_lower_inplace(b, n);
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      double sum = 0.0;
      for (int k = 0; k < n; k++) sum += b[k][i] * b[k][j];
      out[i][j] = sum;
    }
  }
  for (int i = 0; i < n; i++) for (int j = 0; j < i; j++) out[i][j] = out[j][i];
}

static void matvec(const dmat &A, const dvec &x, dvec &z, int n) {
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) sum += A[i][j] * x[j];
    z[i] = sum;
  }
}

static inline double pnorm_c(double y, double m, double s) {   // mombf pnormC clamps
  double z = (y - m) / s;
  if (z < -20.0) return 2.753624e-89;
  if (z > 20.0) return 1 - 2.753624e-89;
  return R::pnorm(y, m, s, 1, 0);
}

static inline double qnorm_c(double cdf, double m, double s) { // mombf qnormC clamps
  if (cdf <= 2.753624e-89) return -20.0 * s + m;
  if (cdf >= 0.99999999999999989) return 8.209536 * s + m;
  return R::qnorm(cdf, m, s, 1, 0);
}

// One draw from N(m, 1) restricted to the union of the disjoint ordered intervals
// [lo[i], hi[i]], i = 0..nint-1 (mombf rnorm_truncMult with n = 1, s = 1)
static double rnorm_trunc_union(const dvec &lo, const dvec &hi, int nint, double m) {
  dvec p0(nint), p1(nint), cump(nint + 1);
  cump[0] = 0.0;
  for (int i = 0; i < nint; i++) {
    p0[i] = pnorm_c(lo[i], m, 1.0);
    p1[i] = pnorm_c(hi[i], m, 1.0);
    cump[i + 1] = std::max(cump[i] + 1.0e-30, cump[i] + p1[i] - p0[i]);
  }
  double u = unif_rand() * cump[nint];
  int j = 0;
  while ((u > cump[j + 1]) && (j < nint - 1)) j++;
  return qnorm_c(p0[j] + u - cump[j], m, 1.0);
}

// One Gibbs sweep for z ~ N(alpha, I) with (D z)_k outside (lower[k], upper[k]) for all k.
// Dz holds D z on input and is kept up to date (mombf rtmvnormOutside_Gibbs).
static void rtmvnorm_outside_gibbs(dvec &z, dvec &Dz, const dvec &alpha, const dmat &D, int p,
                                   const dvec &lower, const dvec &upper) {
  dvec l(p), u(p), lmod(p + 1), umod(p + 1);
  std::vector<int> o(p);
  for (int j = 0; j < p; j++) {
    for (int k = 0; k < p; k++) Dz[k] -= D[k][j] * z[j];
    int nr = 0;
    for (int k = 0; k < p; k++) {
      if (D[k][j] > 0) {
        l[nr] = (lower[k] - Dz[k]) / D[k][j];
        u[nr] = (upper[k] - Dz[k]) / D[k][j];
        nr++;
      } else if (D[k][j] < 0) {
        u[nr] = (lower[k] - Dz[k]) / D[k][j];
        l[nr] = (upper[k] - Dz[k]) / D[k][j];
        nr++;
      }
    }
    if (nr > 0) {
      // excluded regions (l, u), sorted by l and merged; inclusion regions lie in between
      for (int jj = 0; jj < nr; jj++) o[jj] = jj;
      std::sort(o.begin(), o.begin() + nr, [&](int a, int b) { return l[a] < l[b]; });
      int k = 1;                 // number of inclusion intervals - 1 (0-based index of last)
      lmod[0] = l[o[0]]; umod[1] = u[o[0]];
      for (int jj = 1; jj < nr; jj++) {
        if (u[o[jj]] > umod[k]) {
          if (l[o[jj]] <= umod[k]) {
            umod[k] = u[o[jj]];
          } else {
            lmod[k] = l[o[jj]];
            k++;
            umod[k] = u[o[jj]];
          }
        }
      }
      umod[0] = -1.0e20; lmod[k] = 1.0e20;
      // inclusion intervals: [umod[i], lmod[i]], i = 0..k
      z[j] = rnorm_trunc_union(umod, lmod, k + 1, alpha[j]);
      for (int kk = 0; kk < p; kk++) Dz[kk] += z[j] * D[kk][j];
    } else {
      z[j] = alpha[j] + norm_rand();
    }
  }
}

static inline double pen_imom_log(double th, double phi, double tau) {
  // log iMOM density minus log N(0, tau*phi) density (mombf pen_imom)
  double y2 = th * th, s = std::sqrt(tau * phi);
  double dimom = .5 * (std::log(tau) + std::log(phi)) - .5 * kLogPi - std::log(y2) - tau * phi / y2;
  double dnorm = -std::log(2.5066282746310002416123552) - std::log(s) - 0.5 * th * th / (s * s);
  return dimom - dnorm;
}

// Inverse of the iMOM log-penalty: |th| s.t. pen_imom_log(th) = loglambda (mombf invpen_imom_sandwich)
static double invpen_imom(double loglambda, double phi, double tau) {
  const int maxiter = 50; const double ftol = 1.0e-5;
  double b = .5 * (std::log(tau * tau) + 2.0 * std::log(phi) + std::log(2.0)) - loglambda;
  double d = std::sqrt(b * b + 2.0);
  double zcur = tau * phi * (-b + d), thcur = std::sqrt(zcur), fcur = pen_imom_log(thcur, phi, tau);
  double zlow, flow, zup, fup, thlow, thup;
  if (fcur >= loglambda) {
    zlow = .8 * .8 * zcur; thlow = std::sqrt(zlow); flow = pen_imom_log(thlow, phi, tau);
    while (flow >= loglambda) {
      zcur = zlow; thcur = thlow; fcur = flow;
      zlow = .8 * .8 * zcur; thlow = std::sqrt(zlow); flow = pen_imom_log(thlow, phi, tau);
    }
    zup = zcur; fup = fcur;
  } else {
    zup = 1.2 * 1.2 * zcur; thup = std::sqrt(zup); fup = pen_imom_log(thup, phi, tau);
    while (fup <= loglambda) {
      zcur = zup; thcur = thup; fcur = fup;
      zup = 1.2 * 1.2 * zcur; thup = std::sqrt(zup); fup = pen_imom_log(thup, phi, tau);
    }
    zlow = zcur; flow = fcur;
  }
  double err = fcur - loglambda;
  for (int i = 1; i < maxiter && std::fabs(err) > ftol; i++) {
    double slope = (fup - flow) / (zup - zlow);
    zcur = zlow + (loglambda - flow) / slope;
    thcur = std::sqrt(zcur);
    fcur = pen_imom_log(thcur, phi, tau);
    err = fcur - loglambda;
    if (err > 0) { zup = zcur; fup = fcur; } else { zlow = zcur; flow = fcur; }
  }
  return thcur;
}

// Returns the last draw (length ncol(x)); niter Gibbs iterations, the first `burnin` discarded.
// [[Rcpp::export]]
NumericVector bisam_rnlp(NumericVector y, NumericMatrix x, int prior, double tau, double phi,
                         int niter, int burnin, NumericVector thinit) {
  const int n = y.size(), p = x.ncol();
  if (x.nrow() != n) stop("bisam_rnlp: dimensions of y and x do not match");
  if (thinit.size() != p) stop("bisam_rnlp: thinit must have length ncol(x)");
  if (niter - burnin < 1) stop("bisam_rnlp: niter must exceed burnin");
  if (prior != 0 && prior != 1) stop("bisam_rnlp: prior must be 0 (MOM) or 1 (iMOM)");
  const double *X = x.begin(), *Y = y.begin();

  // S = X'X + I/tau, Sinv, chol(Sinv) and its inverse; posterior mean m = Sinv X'y
  dmat S(p, dvec(p)), Sinv(p, dvec(p)), cholSinv(p, dvec(p)), inv_cholSinv(p, dvec(p));
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      double sum = 0.0;
      for (int k = 0; k < n; k++) sum += X[k + (size_t) i * n] * X[k + (size_t) j * n];
      S[i][j] = sum;
    }
  }
  for (int i = 0; i < p; i++) S[i][i] += 1.0 / tau;
  inv_posdef(S, p, Sinv);
  cholesky_lower(Sinv, p, cholSinv);
  cholesky_lower(Sinv, p, inv_cholSinv);
  invert_lower_inplace(inv_cholSinv, p);
  dvec Xty(p), m(p), mortho(p);
  for (int j = 0; j < p; j++) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += X[i + (size_t) j * n] * Y[i];
    Xty[j] = sum;
  }
  matvec(Sinv, Xty, m, p);
  matvec(inv_cholSinv, m, mortho, p);

  // Fixed phi: whitening matrices do not change over iterations
  const double sqrtphi = std::sqrt(phi);
  dmat D(p, dvec(p, 0.0)), Kmat(p, dvec(p, 0.0));
  dvec alpha(p);
  for (int j = 0; j < p; j++) {
    alpha[j] = mortho[j] / sqrtphi;
    for (int k = 0; k <= j; k++) {
      D[j][k] = cholSinv[j][k] * sqrtphi;
      Kmat[j][k] = inv_cholSinv[j][k] / sqrtphi;
    }
  }

  dvec th(thinit.begin(), thinit.end()), lower(p), upper(p), z(p);
  for (int it = 1; it <= niter; it++) {
    for (int i = 0; i < p; i++) {
      if (prior == 0) {           // MOM, r = 1: slice under (th^2 / (phi tau))
        double upperb = th[i] * th[i] / (phi * tau);
        double l = unif_rand() * upperb;
        upper[i] = std::sqrt(l * tau * phi);
      } else {                    // iMOM: slice under the log-penalty
        double l = std::log(unif_rand()) + pen_imom_log(th[i], phi, tau);
        upper[i] = invpen_imom(l, phi, tau);
      }
      lower[i] = -upper[i];
    }
    matvec(Kmat, th, z, p);       // z = K th
    rtmvnorm_outside_gibbs(z, th, alpha, D, p, lower, upper);  // th holds D z
    matvec(D, z, th, p);          // th = D z
  }
  return NumericVector(th.begin(), th.end());
}
