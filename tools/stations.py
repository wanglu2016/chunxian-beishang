import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from onset import onset
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work')
D0 = np.datetime64('1951-01-01')

z = np.load(os.path.join(BASE, 'tmat.npz'))
ids, M = list(z['ids']), z['M'].astype(np.float32)
n, NDAY = M.shape

def fill_gaps(M, maxgap=3):
    M = M.copy()
    for i in range(M.shape[0]):
        x = M[i]
        nanm = np.isnan(x)
        if not nanm.any() or nanm.all():
            continue
        # runs of NaN
        d = np.diff(np.concatenate([[0], nanm.astype(int), [0]]))
        st, en = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
        for a, b in zip(st, en):
            if b - a <= maxgap and a > 0 and b < len(x):
                x[a:b] = np.interp(np.arange(a, b), [a - 1, b], [x[a - 1], x[b]])
    return M

Mf = fill_gaps(M)
np.save(os.path.join(BASE, 'tmat_filled.npy'), Mf)

def di(s):
    return int((np.datetime64(s) - D0).astype(int))

def window(Y):
    a = di(f'{Y-1}-07-28')  # 4 lead days so S is defined from Aug 1
    b = di(f'{Y}-08-31')
    return a, b + 1

def cold_idx(Y):
    a, _ = window(Y)
    return di(f'{Y-1}-11-01') - a, di(f'{Y}-03-31') - a

res = {}
years = list(range(1952, 2021))
for Y in years:
    a, b = window(Y)
    if a < 0:
        a0 = 0
        t = np.concatenate([np.full((n, -a), np.nan, np.float32), Mf[:, 0:b]], axis=1)
    else:
        t = Mf[:, a:b]
    c0, c1 = cold_idx(Y)
    o, f = onset(t, cold=(c0, c1))
    jan1 = di(f'{Y}-01-01') - a
    doy = np.where(f <= 1, o - jan1 + 1, -9999)  # 1 = Jan 1
    res[Y] = (doy, f)

# fallback for Y where previous Aug-Dec is missing (e.g. 1973 after the 1965-72 gap): search from Jan 1
for Y in years:
    doy, f = res[Y]
    a, b = window(Y)
    jan = di(f'{Y}-01-01')
    prev_cov = (~np.isnan(Mf[:, max(a, 0):jan])).sum(axis=1)
    need = (f == 3) & (prev_cov < 30)
    if need.any():
        t = Mf[:, jan - 4:b]
        c1 = di(f'{Y}-03-31') - (jan - 4)
        o, f2 = onset(t, cold=(4, c1))
        # only accept normal onsets (cannot judge 无冬 without the previous autumn)
        okk = need & (f2 == 0)
        d2 = o - 4 + 1
        doy = np.where(okk, d2, doy)
        f = np.where(okk, 0, f)
        res[Y] = (doy, f)

# climatological (常年-style) series: mean annual cycle without Feb 29, from Jul 28 of Y-1 to Aug 31 of Y
def clim(Y0, Y1, minyears):
    stack = []
    for Y in range(Y0, Y1 + 1):
        a, b = window(Y)
        if a < 0:
            continue
        dates = np.arange(np.datetime64(f'{Y-1}-07-28'), np.datetime64(f'{Y}-09-01'))
        keep = ~((dates.astype('datetime64[M]').astype(int) % 12 == 1) & ((dates - dates.astype('datetime64[M]')).astype(int) == 28))
        stack.append(Mf[:, a:b][:, keep])
    X = np.stack(stack)  # years x n x L
    cnt = (~np.isnan(X)).sum(axis=0)
    mean = np.nanmean(X, axis=0)
    mean[cnt < minyears] = np.nan
    return mean

ref_dates = np.arange(np.datetime64('2001-07-28'), np.datetime64('2002-09-01'))  # non-leap reference
def clim_onset(mean):
    c0 = int((np.datetime64('2001-11-01') - ref_dates[0]).astype(int))
    c1 = int((np.datetime64('2002-03-31') - ref_dates[0]).astype(int))
    o, f = onset(mean, cold=(c0, c1))
    jan1 = int((np.datetime64('2002-01-01') - ref_dates[0]).astype(int))
    return np.where(f <= 1, o - jan1 + 1, -9999), f

import warnings; warnings.filterwarnings('ignore')
C = {}
for key, (y0, y1, mny) in {'n9120': (1991, 2020, 20), 'e1': (1957, 1964, 6), 'e2': (2011, 2020, 8), 'n8110': (1981, 2010, 20)}.items():
    C[key] = clim_onset(clim(y0, y1, mny))

out = {'ids': ids, 'years': years,
       'doy': np.array([res[Y][0] for Y in years]).T.tolist(),
       'flag': np.array([res[Y][1] for Y in years]).T.tolist(),
       'clim': {k: [v[0].tolist(), v[1].tolist()] for k, v in C.items()}}
json.dump(out, open(os.path.join(BASE, 'station_onsets.json'), 'w'))
print('ok')
