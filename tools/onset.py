"""Spring-onset (入春) per GB/T 42074-2022 / QX/T 152-2012 logic, vectorised over rows.

Rules implemented
- sliding mean S_j = mean(t_{j-4..j})  (当天及前4天)
- winter process: 5 consecutive S < 10
- spring process: 5 consecutive S >= 10 after the winter process; onset = first day with t >= 10
  among the 9 daily values behind those 5 sliding values (t_{j-4..j+4})
- secondary judgement: if a later 5-day S<10 process occurs, count S>=10 vs S<10 days between the
  candidate onset and the next spring process; if satisfied < unsatisfied, move to the next process
  (repeat)
- no winter process in the cold half-year -> 无冬; onset proxy = day of minimum S (Nov 1 - Mar 31)
- no spring process before window end -> 无春
"""
import numpy as np

TH = 10.0

def sliding(t):
    """t: (n, L) float with NaN; trailing 5-day mean, NaN if any NaN."""
    n, L = t.shape
    S = np.full((n, L), np.nan, dtype=np.float32)
    c = np.nancumsum(np.nan_to_num(t, nan=0.0), axis=1, dtype=np.float64)
    nn = np.cumsum(np.isnan(t), axis=1)
    c = np.concatenate([np.zeros((n, 1)), c], axis=1)
    nn = np.concatenate([np.zeros((n, 1), dtype=int), nn], axis=1)
    s5 = (c[:, 5:] - c[:, :-5]) / 5.0
    bad = (nn[:, 5:] - nn[:, :-5]) > 0
    s5[bad] = np.nan
    S[:, 4:] = s5
    return S

def run5(b):
    """b: (n, L) bool -> r[:, k] = b[k..k+4] all True (False near end)."""
    n, L = b.shape
    r = np.zeros_like(b)
    r[:, :L - 4] = b[:, :L - 4] & b[:, 1:L - 3] & b[:, 2:L - 2] & b[:, 3:L - 1] & b[:, 4:]
    return r

def first_true_after(r, start):
    """r: (n, L) bool, start: (n,) int -> first index >= start where r, else -1."""
    n, L = r.shape
    idx = np.arange(L)[None, :]
    m = r & (idx >= start[:, None])
    has = m.any(axis=1)
    out = np.where(has, m.argmax(axis=1), -1)
    return out

def onset(t, cold=(92, 243), strict_nan=True):
    """
    t: (n, L) daily mean temps for a window starting Aug 1 of previous year (index 0).
    cold: index range of the cold half-year (default Nov 1 .. Mar 31 for a window from Aug 1).
    Returns onset index (into window), flag: 0 normal, 1 无冬, 2 无春, 3 NA.
    """
    n, L = t.shape
    S = sliding(t)
    isn = np.isnan(S)
    isn[:, :4] = False  # first 4 have no sliding value by construction (handled as unknown below)
    A = S >= TH
    B = S < TH
    A5, B5 = run5(A), run5(B)
    ar = np.arange(L)[None, :]
    nanc = np.concatenate([np.zeros((n, 1), int), np.cumsum(np.isnan(S[:, :]) & (ar >= 4), axis=1)], axis=1)
    def nan_between(a, b):  # count NaN sliding values in [a, b]
        a = np.clip(a, 0, L - 1); b = np.clip(b, 0, L - 1)
        return nanc[np.arange(n), b + 1] - nanc[np.arange(n), a]

    out = np.full(n, -9999, dtype=np.int32)
    flag = np.full(n, 3, dtype=np.int8)
    kw = first_true_after(B5, np.full(n, 4))
    # secondary judgement for the winter process too: a cold snap followed by a longer warm
    # spell is not the start of winter (otherwise a November cold snap in the south would
    # make a December "spring")
    csA0 = np.concatenate([np.zeros((n, 1), int), np.cumsum(A, axis=1)], axis=1)
    csB0 = np.concatenate([np.zeros((n, 1), int), np.cumsum(B, axis=1)], axis=1)
    rows0 = np.arange(n)
    act = kw >= 0
    for _ in range(8):
        m1 = first_true_after(A5, np.where(act, kw + 5, L))
        k2 = first_true_after(B5, np.where(act & (m1 >= 0), m1 + 5, L))
        res = act & (m1 >= 0) & (k2 >= 0) & (k2 <= cold[1])
        ncold = csB0[rows0, np.clip(k2, 0, L)] - csB0[rows0, np.clip(kw, 0, L)]
        nwarm = csA0[rows0, np.clip(k2, 0, L)] - csA0[rows0, np.clip(kw, 0, L)]
        mv = res & (ncold < nwarm)
        if not mv.any():
            break
        kw = np.where(mv, k2, kw)
        act = mv
    # winter must start before the end of the cold half-year to count as this winter
    kw = np.where(kw > cold[1], -1, kw)
    # no winter
    now = kw < 0
    cold_nan = nan_between(np.full(n, cold[0]), np.full(n, cold[1]))
    nw_ok = now & (cold_nan == 0)
    if nw_ok.any():
        Sc = S[nw_ok, cold[0]:cold[1] + 1]
        out[nw_ok] = cold[0] + np.nanargmin(Sc, axis=1)
        flag[nw_ok] = 1
    # winter present: find spring process
    w = ~now
    j = np.where(w, first_true_after(A5, np.maximum(kw + 5, 4)), -1)
    # NaN guard between winter start and spring process end
    nb = nan_between(np.maximum(kw, 4), np.where(j >= 0, j + 4, L - 1))
    ok = w & (nb == 0)
    nos = ok & (j < 0)
    flag[nos] = 2
    cand = ok & (j >= 0)

    def onset_day(jj):
        # first t>=10 within t[jj-4 .. jj+4]
        rows = np.arange(n)
        best = np.full(n, -1)
        for off in range(8, -1, -1):
            d = np.clip(jj - 4 + off, 0, L - 1)
            hit = t[rows, d] >= TH
            best = np.where(hit, d, best)
        return np.where(best >= 0, best, jj)

    jc = np.where(cand, j, 0)
    dc = onset_day(jc)
    csA = np.concatenate([np.zeros((n, 1), int), np.cumsum(A, axis=1)], axis=1)
    csB = np.concatenate([np.zeros((n, 1), int), np.cumsum(B, axis=1)], axis=1)
    active = cand.copy()
    for _ in range(6):
        k = first_true_after(B5, jc + 1)
        k = np.where(active, k, -1)
        has_dip = active & (k >= 0)
        m = first_true_after(A5, np.where(has_dip, k + 5, L))
        m = np.where(has_dip, m, -1)
        resolve = has_dip & (m >= 0)
        rows = np.arange(n)
        nsat = csA[rows, np.clip(m, 0, L)] - csA[rows, dc]
        nun = csB[rows, np.clip(m, 0, L)] - csB[rows, dc]
        move = resolve & (nsat < nun)
        if not move.any():
            break
        jc = np.where(move, m, jc)
        dc = np.where(move, onset_day(jc), dc)
        active = move
    # the 30 days after the accepted spring process must be observed, otherwise the
    # secondary judgement could not have been made -> NA
    # (and, for early candidates, the rest of the cold half-year too)
    seen = nan_between(jc, np.minimum(np.maximum(jc + 34, cold[1]), L - 1)) == 0
    acc = cand & seen
    out[acc] = dc[acc]
    flag[acc] = 0
    return out, flag
