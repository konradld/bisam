"""
verify_appendix_S1.py -- re-runnable checks for every numerical claim in Appendix S1 (v3).

  python3 verify_appendix_S1.py          (needs numpy, scipy, mpmath; ~3 min)

Each block prints PASS/FAIL and the numbers quoted in the appendix.
"""
import numpy as np, mpmath as mp
from scipy.special import gamma as Gam
from numpy.polynomial.hermite_e import hermegauss
mp.mp.dps = 50
rng = np.random.default_rng(0)
OK = lambda b: "PASS" if b else "**FAIL**"

def u(s, T): return (np.arange(1, T+1) >= s).astype(float)
def V(m, T): return np.column_stack([np.ones(T)]+[u(s, T) for s in sorted(m)])
def Pr(A):
    Q, _ = np.linalg.qr(A); return Q @ Q.T

# ---------------------------------------------------------------- 1. Prop S1.1
err = 0.
for _ in range(300):
    T = int(rng.integers(12, 60))
    m = sorted(rng.choice(np.arange(3, T), size=int(rng.integers(0, 6)), replace=False))
    s = int(rng.integers(3, T))
    if s in m: continue
    mm = [x for x in m if x != s]
    H = Pr(V(mm, T)) if mm else Pr(np.ones((T, 1)))
    ut = (np.eye(T)-H) @ u(s, T)
    a = max([x for x in mm if x < s]+[1]); b = min([x for x in mm if x > s]+[T+1])
    ind = lambda lo, hi: ((np.arange(1, T+1) >= lo) & (np.arange(1, T+1) < hi)).astype(float)
    ns = (s-a)*(b-s)/(b-a)
    err = max(err, np.abs(ut-(ind(s, b)-(b-s)/(b-a)*ind(a, b))).max(), abs(ut @ ut-ns))
    err = max(err, max(0, .5*min(s-a, b-s)-ns), max(0, ns-min(s-a, b-s)))
    for sp in range(3, T):
        cl = 0. if not (a < sp < b) else (min(s, sp)-a)*(b-max(s, sp))/(b-a)
        err = max(err, abs(ut @ u(sp, T)-cl))
print(f"1. Prop S1.1 (i)-(iii) and the n_s sandwich: max error {err:.2e}  {OK(err<1e-10)}")

# ------------------------------------------------- 2. Lemma S1.3 energy bound
def energy(m, mu, T): return float(np.sum(((np.eye(T)-Pr(V(m, T))) @ mu)**2))
worst = np.inf
for _ in range(400):
    T = int(rng.integers(80, 160)); k0 = int(rng.integers(1, 4)); R = int(rng.integers(2, 8))
    d = None
    for _ in range(200):
        c = np.sort(rng.choice(np.arange(3, T-1), size=k0, replace=False))
        if np.diff(np.concatenate(([1], c, [T+1]))).min() >= 2*R: d = c; break
    if d is None: continue
    g = np.array([(-1)**i for i in range(k0)])*rng.uniform(.8, 3., size=k0)
    mu = sum(gi*u(si, T) for gi, si in zip(g, d))
    m = sorted(set(rng.choice(np.arange(3, T), size=int(rng.integers(0, 7)), replace=False).tolist()))
    U = sum(1 for si in d if min([abs(si-x) for x in m], default=T) > R)
    if U == 0: continue
    worst = min(worst, energy(m, mu, T)/(.5*np.abs(g).min()**2*R*U))
e_bd = energy([60], u(4, 100), 100)
print(f"2. Lemma S1.3 with boundary-inclusive Delta_T: min ratio {worst:.3f}  {OK(worst>=1)}")
print(f"   without it (T=100, s0=4, m={{60}}, R=12): energy {e_bd:.2f} vs 6.00 -> bound fails  "
      f"{OK(e_bd < 6)}")

# ------------------------------------------- 3. Lemma S1.4 exact local factor
err = 0.
for _ in range(200):
    T = int(rng.integers(30, 90)); g = float(rng.choice([1., 10., float(T)])); sig = 1.
    m = sorted(set(rng.choice(np.arange(3, T), size=int(rng.integers(1, 4)), replace=False).tolist()))
    s = int(rng.integers(3, T))
    if s in m: continue
    z = rng.normal(size=T)*1.3+2*u(T//2, T)
    M1 = np.eye(T)-Pr(np.ones((T, 1))); zc = M1 @ z
    f = lambda mm: (-.5*len(mm)*np.log(1+g)
                    - (np.sum(zc**2)-g/(1+g)*np.sum((Pr(np.column_stack([M1 @ u(x, T) for x in mm])) @ zc)**2))/(2*sig**2)
                    if mm else -np.sum(zc**2)/(2*sig**2))
    ut = (np.eye(T)-Pr(V(m, T))) @ u(s, T); t2 = (ut @ z)**2/(ut @ ut)/sig**2
    err = max(err, abs(f(sorted(m+[s]))-f(m)-(-.5*np.log(1+g)+g/(1+g)*t2/2)))
print(f"3. Lemma S1.4 (local factor exactly scalar in t_s): max error {err:.2e}  {OK(err<1e-9)}")

# ------------------------------- 4. Lemma S1.5: beta_k, side condition, bounds
bstat = lambda k: (8*k)**(-k/(k+1))+(8*k)**(1/(k+1))/8
bshrp = lambda k: (2*k)**(-k/(k+1))+(2*k)**(1/(k+1))/2
ok = True
for k in [.5, 1, 2, 3]:
    w = np.linspace(1e-4, 5, 500000)
    for lam in [10, 1000]:
        ok &= abs((w**(-2*k)+lam*w**2/8).min()-bstat(k)*lam**(k/(k+1))) < 1e-3
print(f"4a. beta_k = min_w{{|w|^-2k + lam w^2/8}} lam^-k/(k+1): {OK(ok)}")
sc = {k: ((2*bstat(k))**(-k) >= bstat(k)/2) for k in [.25, .5, 1, 2, 3, 5]}
print(f"4b. side condition (2 beta_k)^-k >= beta_k/2 for k in [1/4,5]: {OK(all(sc.values()))}")
print(f"    beta_1={bstat(1):.5f}, beta*_1={bshrp(1):.5f} (=sqrt2 {OK(abs(bshrp(1)-2**.5)<1e-12)}),"
      f" conservatism factor {bshrp(1)/(.5*bstat(1)):.2f}")

def BF(t, lam, k=1, nu=1):
    t, lam, k, nu = map(mp.mpf, (t, lam, k, nu)); c = k/mp.gamma(nu/(2*k))
    f = lambda w: c*mp.e**(-(nu+1)*mp.log(abs(w))-abs(w)**(-2*k)+t*w*mp.sqrt(lam)-lam*w**2/2)
    return (mp.quad(f, [mp.mpf('1e-8'), .05, .5, 2, 20, mp.inf])
            + mp.quad(f, [-mp.inf, -20, -2, -.5, -.05, mp.mpf('-1e-8')]))
Amax = 0.
for k, nu in [(1, 1), (1, 3), (2, 1)]:
    for lam in [1, 10, 100, 2000]:
        t0 = mp.sqrt(bstat(k)*mp.mpf(lam)**(k/(k+1))/2)
        for fr in [0, .5, 1.]:
            r = (max(BF(t0*fr, lam, k, nu), BF(-t0*fr, lam, k, nu))
                 / (mp.mpf(lam)**(nu/2)*mp.e**((t0*fr)**2/2-bstat(k)*mp.mpf(lam)**(k/(k+1))/2)))
            Amax = max(Amax, float(r))
print(f"4c. Lemma S1.5(i) holds with A={Amax:.3g} <= 1: {OK(Amax<=1)}")
slack = 1e9
for k, nu in [(1, 1), (2, 3)]:
    for lam in [1, 50, 5000]:
        for t in [-4, 0.3, 3, 20]:
            wst = t/np.sqrt(lam)
            gr = [w for w in np.linspace(wst-lam**-.5, wst+lam**-.5, 41) if abs(w) > 1e-12]
            mn = min(float(mp.log(k/mp.gamma(mp.mpf(nu)/(2*k)))-(nu+1)*mp.log(abs(w))-abs(w)**(-2*k)) for w in gr)
            slack = min(slack, float(mp.log(BF(t, lam, k, nu)))-(t**2/2-.5*np.log(lam)-.5+np.log(2)+mn))
print(f"4d. Lemma S1.5(ii) slack {slack:.3f} >= 0: {OK(slack>=0)}")
r = float(mp.log(BF(0, 10000, 1, 1)))/100
print(f"4e. sharp null rate log BF/lam^(1/2) at lam=1e4: {r:.4f} vs -beta*_1={-bshrp(1):.4f}  {OK(abs(r+bshrp(1))<.02)}")

# ----------------------------- 5. Prop S1.6 Step 1: the exponent is an identity
err = 0.
for _ in range(200):
    T = int(rng.integers(30, 80)); sig = float(rng.uniform(.5, 2))
    d = sorted(rng.choice(np.arange(5, T-3), size=2, replace=False)); g0 = rng.normal(size=2)*1.5
    mu = sum(gi*u(si, T) for gi, si in zip(g0, d)); eps = rng.normal(scale=sig, size=T); z = mu+eps
    m = [x for x in sorted(set(rng.choice(np.arange(3, T), size=3, replace=False).tolist())) if x not in d]
    M = np.eye(T)-Pr(V(m, T)); XJ = np.column_stack([M @ u(si, T) for si in d])
    eta = rng.normal(size=2)*.05; r = M @ z
    lhs = (r @ (XJ @ (g0+eta)))/sig**2-np.sum((XJ @ (g0+eta))**2)/(2*sig**2)
    U = float(np.sum((M @ mu)**2))
    rhs = (U/(2*sig**2)+((M @ mu) @ eps)/sig**2+((M @ eps) @ (XJ @ eta))/sig**2
           - np.sum((XJ @ eta)**2)/(2*sig**2))
    err = max(err, abs(lhs-rhs))
print(f"5. Prop S1.6 Step 1 exponent decomposition: max error {err:.2e}  {OK(err<1e-9)}")

# ------------------------------------------------- 6. Corollary S1.10 arithmetic
t0sq = lambda N, T, h: (np.sqrt(2*np.log(2*N*T**3/0.05))+2*np.sqrt(h))**2
print("6. Corollary S1.10: log10 of the RHS of (S1.gaugebound), tau=3.32, n_low=T/4, h=log(NT)")
print(f"   {'T':>7} {'N':>3} {'iMOM k=1':>10} {'pMOM k=1':>10} {'pMOM k=3':>10} {'local':>8}")
for T in [1e4, 1e5]:
    for N in [1, 30]:
        nl = T/4; e = t0sq(N, T, np.log(N*T))/2/np.log(10); lam = 3.32*nl
        print(f"   {T:>7.0e} {N:>3} {np.log10(lam**.5)+e-.5*bstat(1)*lam**.5/np.log(10):>10.1f}"
              f" {e-1.5*np.log10(nl):>10.1f} {e-3.5*np.log10(nl):>10.1f} {e-.5*np.log10(nl):>8.1f}")
print("   (v2's accounting omits e^{t0^2/2}; then pMOM k=1 reads -5.1 at T=1e4, N=1)")

# --------------------------------------------------- 7. window rule and the cap
bad = []
for kap in [1/3, 1/2, .1, .01]:
    tm = int(np.ceil(1+1/kap)-1)
    if not (kap < (1/(tm-1) if tm > 1 else np.inf)): bad.append(kap)
print(f"7. cap t_max = ceil(1+c0/kappa)-1 implies kappa < c0/(t_max-1): {OK(not bad)}")
print(f"   defaults (1,1,1/3): t_max=3, P^1={0.5:.3f}, P^2={(1+1/3)/2:.3f}, P^3={(1+2/3)/2:.3f}")
for N, T in [(15, 30), (30, 60)]:
    print(f"   potency needs t_max >> log K = {np.log(N*(T-3)):.1f} (N={N},T={T}) but default is 3")

# --------------------- 8. clustering: n_s for adjacent dates, and lambda_min
T = 500
ns = [float(((np.eye(T)-Pr(V([50, 250+gp], T))) @ u(250, T)) @ ((np.eye(T)-Pr(V([50, 250+gp], T))) @ u(250, T)))
      for gp in (1, 2, 3)]
print(f"8. n_s at separations 1,2,3: {ns[0]:.3f}, {ns[1]:.3f}, {ns[2]:.3f}  {OK(ns[0]<1)}")
Xm = np.column_stack([u(s, 400) for s in np.linspace(3, 399, 40).astype(int)])
lm = np.linalg.eigvalsh(Xm.T @ Xm).min()
print(f"   lambda_min(X_m'X_m)={lm:.3f} <= T/(2(|m|-1))={400/78:.2f}  {OK(lm<=400/78)}")

# ------------------ 9. Prop S1.14(i): displacement energy uses the segment length
T, s0, a, b = 400, 200, 100, 300
for r in (1, 3, 10):
    ex = energy([a, s0+r, b], u(s0, T), T); n1 = s0+r-a
    print(f"9. r={r:>2}: exact {ex:.4f}, segment formula {r*(n1-r)/n1:.4f}, v2's (b-a) formula "
          f"{r*(b-a-r)/(b-a):.4f}  {OK(abs(ex-r*(n1-r)/n1)<1e-6)}")