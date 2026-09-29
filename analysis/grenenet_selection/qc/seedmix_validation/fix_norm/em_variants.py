"""EM variants for the founder-completeness-bias fix (imported, self-contained).

Do NOT edit src/kmate. This is a standalone copy of the M-step with a
per-founder normalization toggle.

Model (noiseless): c_k = sum_f h_f K[f,k],  mu_k = sum_f h_f K[f,k].

E-step responsibilities:  E[n_fk] = c_k * h_f K[f,k] / mu_k.

M-step variants:
  'multinomial' (PRODUCTION current):  h_new_f = (sum_k E[n_fk]) / total_c
        -> h_new_f propto h_f * sum_k K[f,k] c_k/mu_k
        -> at h=h_true this is h_f * Kf_weighted (NOT a fixed point):
           over-credits founders with MORE k-mers (assembly completeness).
  'poisson'  (THE FIX):  h_new_f = (sum_k E[n_fk]) / (sum_k omega_k K[f,k])
        -> divide each founder's evidence by its own (weighted) k-mer content.
        -> at h=h_true, c=mu so numerator = h_f * sum_k omega_k K[f,k] = h_f*Kf_w,
           denominator = Kf_w  => h_new_f = h_f. EXACT fixed point for ANY h_true
           and ANY omega. Completeness cancels. Renormalize to simplex.

omega_k per-k-mer weights (graded down-weighting, optional):
  None       -> omega=1
  '1/ac'     -> 1/carriers
  '1/sqrtac' -> 1/sqrt(carriers)
  '1/mb'     -> 1/bubble-size (production bubble de-replication)
"""
import numpy as np


def make_omega(kind, ac, bubble_id=None):
    if kind is None:
        return None
    if kind == "1/ac":
        return (1.0 / np.maximum(ac, 1)).astype(np.float32)
    if kind == "1/sqrtac":
        return (1.0 / np.sqrt(np.maximum(ac, 1))).astype(np.float32)
    if kind == "1/mb":
        assert bubble_id is not None
        # bubble size = #kmers sharing a bubble_id
        uniq, inv, cnt = np.unique(bubble_id, return_inverse=True, return_counts=True)
        return (1.0 / cnt[inv]).astype(np.float32)
    raise ValueError(kind)


def solve_em(counts, K_csr, mode="multinomial", omega=None,
             h_init=None, max_iter=300, tol=1e-9, verbose=False):
    """K_csr: F x Kn scipy csr float32. counts: Kn vector."""
    F, Kn = K_csr.shape
    counts = counts.astype(np.float32)
    wc = counts if omega is None else (omega.astype(np.float32) * counts)
    total_c = float(wc.sum())

    # per-founder (weighted) k-mer content Kf_w = sum_k omega_k K[f,k]
    if omega is None:
        Kf_w = np.asarray(K_csr.sum(axis=1)).ravel().astype(np.float32)
    else:
        Kf_w = np.asarray(K_csr @ omega.astype(np.float32)).ravel().astype(np.float32)
    Kf_w = np.maximum(Kf_w, 1e-12)

    h = np.full(F, 1.0 / F, np.float32) if h_init is None else h_init.astype(np.float32).copy()
    Kt = K_csr.T  # csc VIEW (shares arrays, no copy) -> Kt @ h == h @ K

    hist = []
    for it in range(max_iter):
        mu = np.maximum(Kt @ h, np.float32(1e-9))     # Kn
        cw = wc / mu                                   # Kn
        em = h * (K_csr @ cw)                          # F
        if mode == "multinomial":
            h_new = em / max(total_c, 1e-12)
        elif mode == "poisson":
            h_new = em / Kf_w
        else:
            raise ValueError(mode)
        h_new = h_new / h_new.sum()
        d = float(np.linalg.norm(h_new - h))
        hist.append(d)
        h = h_new
        if verbose and it % 25 == 0:
            print(f"    it{it} d={d:.2e} min={h.min():.2e} max={h.max():.4f}")
        if d < tol:
            break
    return h.astype(np.float64), {"iters": it + 1, "conv": d < tol, "last_d": d, "Kf_w": Kf_w}
