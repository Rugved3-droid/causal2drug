#!/usr/bin/env python
"""SuSiE-RSS + coloc.susie -- credible-set (LD-aware) colocalization.

METHOD CHOICE AND WHY
---------------------
Single-variant coloc (coloc.abf) compares one variant at a time and its evidence scales
as ~(1/2)z^2, so when two traits' observed peaks sit at different positions inside one LD
block the log-ABF gap explodes and PP4 collapses -- measured at HMGCR-LDL-C as a gap of
32.23 log units giving PP4 = 1e-11 (deliverables/d1d).

coloc.susie (Wallace 2021) is the standard fix and is what is implemented here. It runs
SuSiE-RSS (Zou et al. 2022) on each trait's z-scores with an LD matrix, which decomposes
the region into L single-effect components with per-SNP posterior inclusion probabilities,
then applies the coloc.abf arithmetic to each PAIR of single effects using their log Bayes
factor vectors. Because a single effect's lBF vector is spread across all variants in LD
with the causal one, it is insensitive to which variant happens to peak -- exactly the
failure mode being fixed.

SuSiE-RSS was preferred over conditional/stepwise coloc because it never inverts the LD
matrix (it only multiplies by it), so it is numerically stable at the n_snps > n_samples
ratios this project has (271 SNPs from 525 EUR samples), where an LD-matrix inversion
would be the dominant source of error.

Model (summary-statistic SuSiE): z ~ N(R b, R), fitted by IBSS.
"""
import numpy as np


# ----------------------------- SuSiE-RSS -----------------------------
def _single_effect(rbar, prior_v, pi, s2=1.0):
    """Posterior for one single effect given residual z-scores rbar.

    s2 is the RESIDUAL VARIANCE. Omitting it (fixing s2=1) was a real defect: in a region
    with z up to ~50, the residual after removing a component is itself large, every
    subsequent single effect becomes a delta function, and IBSS carves the region into
    spurious 'signals'. That produced a credible set for the outcome topping at the
    exposure's lead variant with PIP 1.000 while the true GWAS lead had PIP 5e-87 -- an
    artefactual PP4 = 1. Estimating s2 is the standard SuSiE remedy.
    """
    v_over_s2 = prior_v / s2
    lbf = (-0.5 * np.log1p(v_over_s2)
           + 0.5 * (rbar ** 2 / s2) * v_over_s2 / (1.0 + v_over_s2))
    w = lbf + np.log(pi)
    w -= w.max()
    alpha = np.exp(w)
    alpha /= alpha.sum()
    post_var = prior_v / (1.0 + v_over_s2)
    post_mean = post_var * rbar / s2
    return alpha, post_mean, lbf


def _lse(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


PRIOR_GRID = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0])


def _estimate_prior_v(rbar, pi, s2=1.0, grid=PRIOR_GRID):
    """Empirical-Bayes prior variance: pick the grid value maximising the single-effect
    marginal likelihood. Removes the prior variance as a hand-set free parameter --
    important because the effects here are on the z-score scale, where a small prior
    (the value used in a first draft) makes every PIP diffuse and suppresses credible
    sets entirely."""
    best_v, best_ll = 0.0, _lse(np.zeros_like(rbar) + np.log(pi))  # v=0 => null
    for v in grid * s2:
        r = v / s2
        lbf = -0.5 * np.log1p(r) + 0.5 * (rbar ** 2 / s2) * r / (1.0 + r)
        ll = _lse(lbf + np.log(pi))
        if ll > best_ll:
            best_ll, best_v = ll, v
    return float(best_v)


def susie_rss(z, R, L=5, prior_v=None, n_iter=100, tol=1e-6, pi=None):
    """Fit SuSiE-RSS. prior_v=None estimates it by empirical Bayes (default)."""
    z = np.asarray(z, dtype=float)
    p = len(z)
    if pi is None:
        pi = np.full(p, 1.0 / p)
    alpha = np.zeros((L, p))
    mu = np.zeros((L, p))
    lbf = np.zeros((L, p))
    b = np.zeros((L, p))          # alpha * mu, the posterior mean effect of each component
    vs = np.zeros(L)
    s2 = 1.0
    prev = None
    for _ in range(n_iter):
        for l in range(L):
            others = b.sum(axis=0) - b[l]
            rbar = z - R @ others
            v = (_estimate_prior_v(rbar, pi, s2) if prior_v is None else prior_v)
            vs[l] = v
            if v <= 0:            # component explains nothing: leave it null
                alpha[l] = pi; mu[l] = 0.0; lbf[l] = 0.0; b[l] = 0.0
                continue
            a, m, lb = _single_effect(rbar, v, pi, s2)
            alpha[l], mu[l], lbf[l] = a, m, lb
            b[l] = a * m
        fit = R @ b.sum(axis=0)
        resid = z - fit
        # residual variance, floored at 1 because z-scores have unit null variance
        s2 = max(1.0, float(np.mean(resid ** 2)))
        obj = float(np.sum(resid ** 2))
        if prev is not None and abs(prev - obj) < tol:
            break
        prev = obj
    return dict(alpha=alpha, lbf=lbf, b=b, L=L, prior_v=vs, s2=s2)


def credible_sets(fit, R, coverage=0.95, min_abs_corr=0.5):
    """95% credible sets with a purity filter, as in SuSiE."""
    out = []
    for l in range(fit["L"]):
        a = fit["alpha"][l]
        order = np.argsort(-a)
        c = np.cumsum(a[order])
        k = int(np.searchsorted(c, coverage) + 1)
        idx = np.sort(order[:k])
        if len(idx) > 1:
            sub = np.abs(R[np.ix_(idx, idx)])
            purity = sub[np.triu_indices(len(idx), 1)].min()
        else:
            purity = 1.0
        if purity >= min_abs_corr:
            out.append(dict(component=l, idx=idx, size=len(idx), purity=float(purity),
                            top=int(order[0]), top_pip=float(a[order[0]])))
    # drop duplicate components that converged onto the same signal
    seen, uniq = set(), []
    for cs in out:
        key = cs["top"]
        if key in seen:
            continue
        seen.add(key)
        uniq.append(cs)
    return uniq


# --------------------------- coloc.susie ---------------------------
def _ldiff(a, b):
    return a + np.log1p(-np.exp(b - a)) if b < a else -np.inf


def coloc_from_lbf(lbf1, lbf2, p1=1e-4, p2=1e-4, p12=1e-5):
    """coloc.abf arithmetic applied to two log-Bayes-factor vectors."""
    s1, s2 = _lse(lbf1), _lse(lbf2)
    s12 = _lse(lbf1 + lbf2)
    ls = np.array([0.0, np.log(p1) + s1, np.log(p2) + s2,
                   np.log(p1) + np.log(p2) + _ldiff(s1 + s2, s12),
                   np.log(p12) + s12])
    d = _lse(ls)
    pp = np.exp(ls - d)
    return dict(PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3], PP4=pp[4],
                s1=s1, s2=s2, s12=s12, gap=s1 + s2 - s12)


class ColocArtefact(Exception):
    """Raised when the winning credible-set pair carries the over-fitting signature.

    Four PP4 ~ 1 values were produced by this project and withdrawn, every one with the
    same signature: a SuSiE component from each trait converging onto the SAME variant,
    which then colocalizes with itself by construction. Twice the variant was not even a
    lead for either trait (HMGCR/UKB: rs1544755 x rs1544755). Detecting this by eye is not
    a control; it is now an assertion that raises.
    """


def _artefact_reason(pair, lead1, lead2):
    if pair["cs1_top"] == pair["cs2_top"]:
        return ("both credible sets top at the SAME variant (index %d) -- a component "
                "colocalizing with itself" % pair["cs1_top"])
    if pair["cs1_top"] != lead1 and pair["cs2_top"] != lead2:
        return ("winning pair tops (%d, %d) match neither trait's lead (%d, %d)"
                % (pair["cs1_top"], pair["cs2_top"], lead1, lead2))
    return None


def coloc_susie(z1, z2, R, L=5, prior_v=None, strict=True, **pri):
    """Full credible-set colocalization. Returns per-CS-pair results.

    strict=True (default) raises ColocArtefact if the BEST pair carries the over-fitting
    signature. It never silently returns such a number.
    """
    z1 = np.asarray(z1, dtype=float)
    z2 = np.asarray(z2, dtype=float)
    f1 = susie_rss(z1, R, L=L, prior_v=prior_v)
    f2 = susie_rss(z2, R, L=L, prior_v=prior_v)
    cs1 = credible_sets(f1, R)
    cs2 = credible_sets(f2, R)
    lead1 = int(np.argmax(np.abs(z1)))
    lead2 = int(np.argmax(np.abs(z2)))
    pairs = []
    for a in cs1:
        for b in cs2:
            r = coloc_from_lbf(f1["lbf"][a["component"]], f2["lbf"][b["component"]], **pri)
            r.update(cs1_component=a["component"], cs2_component=b["component"],
                     cs1_size=a["size"], cs2_size=b["size"],
                     cs1_purity=a["purity"], cs2_purity=b["purity"],
                     cs1_top=a["top"], cs2_top=b["top"])
            r["artefact_reason"] = _artefact_reason(r, lead1, lead2)
            pairs.append(r)
    best = max(pairs, key=lambda r: r["PP4"]) if pairs else None
    if strict and best is not None and best["artefact_reason"]:
        raise ColocArtefact(
            "PP4=%.4g REJECTED: %s" % (best["PP4"], best["artefact_reason"]))
    return dict(n_cs1=len(cs1), n_cs2=len(cs2), pairs=pairs, best=best,
                cs1=cs1, cs2=cs2, lead1=lead1, lead2=lead2)


def ld_condition_number(R):
    ev = np.linalg.eigvalsh(R)
    lo, hi = float(ev.min()), float(ev.max())
    return dict(min_eigenvalue=lo, max_eigenvalue=hi,
                condition_number=(hi / lo if lo > 1e-12 else float("inf")))
