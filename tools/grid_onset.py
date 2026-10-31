"""Interpolate daily mean temperature to a 20 km grid (daily regression on latitude+elevation,
then inverse-distance weighting of residuals), then apply the same onset rule per grid cell."""
import os, sys, json, numpy as np, warnings
sys.path.insert(0, os.path.dirname(__file__))
from onset import onset
from geo import proj
warnings.filterwarnings('ignore')
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work')
D0 = np.datetime64('1951-01-01')

Mf = np.load(os.path.join(BASE, 'tmat_filled.npy'))
ids = list(np.load(os.path.join(BASE, 'tmat.npz'))['ids'])
meta = json.load(open(os.path.join(BASE, 'meta.json')))
G = np.load(os.path.join(BASE, 'grid.npz'))
near = G['near']
cx, cy = G['GX'][near], G['GY'][near]
clat, cz = G['glat'][near], G['gz'][near]
NC = near.sum()

sidx = np.array([i for i, s in enumerate(ids) if s in meta])
slat = np.array([meta[ids[i]]['lat'] for i in sidx])
slon = np.array([meta[ids[i]]['lon'] for i in sidx])
sz = np.array([meta[ids[i]]['z'] for i in sidx], float)
sx, sy = proj(slon, slat)
K = 10

def di(s):
    return int((np.datetime64(s) - D0).astype(int))

def interp(T):
    """T: (ns, L) station series (rows aligned with sidx). Returns (NC, L) grid series."""
    ns, L = T.shape
    cov = (~np.isnan(T)).mean(axis=1)
    use = cov >= 0.7
    Tu, la, zz, xx, yy = T[use], slat[use], sz[use], sx[use], sy[use]
    # k nearest stations for every cell
    d2 = (cx[:, None] - xx[None, :]) ** 2 + (cy[:, None] - yy[None, :]) ** 2
    order = np.argsort(d2, axis=1)[:, :K]
    W = np.zeros_like(d2, dtype=np.float32)
    rr = np.arange(NC)[:, None]
    W[rr, order] = 1.0 / np.maximum(d2[rr, order], 15.0 ** 2)
    out = np.full((NC, L), np.nan, np.float32)
    X = np.c_[np.ones(len(la)), la, zz / 1000.0]
    Xc = np.c_[np.ones(NC), clat, cz / 1000.0]
    V = ~np.isnan(Tu)
    coefs = np.full((L, 3), np.nan)
    for d in range(L):
        v = V[:, d]
        if v.sum() < 40:
            continue
        coefs[d] = np.linalg.lstsq(X[v], Tu[v, d], rcond=None)[0]
    pred_s = X @ np.nan_to_num(coefs.T)  # ns x L
    res = np.where(V, Tu - pred_s, 0.0).astype(np.float32)
    num = W @ res
    den = W @ V.astype(np.float32)
    ri = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
    out = (Xc @ np.nan_to_num(coefs.T)).astype(np.float32) + ri
    out[:, np.isnan(coefs[:, 0])] = np.nan
    return out, int(use.sum()), coefs

def window(Y):
    return di(f'{Y-1}-07-28'), di(f'{Y}-08-31') + 1

results = {}
def run_year(Y):
    a, b = window(Y)
    T = Mf[sidx, max(a, 0):b]
    if a < 0:
        T = np.concatenate([np.full((len(sidx), -a), np.nan, np.float32), T], axis=1)
    g, nuse, _ = interp(T)
    c0, c1 = di(f'{Y-1}-11-01') - a, di(f'{Y}-03-31') - a
    o, f = onset(g, cold=(c0, c1))
    jan1 = di(f'{Y}-01-01') - a
    doy = np.where(f <= 1, o - jan1 + 1, -9999)
    return doy, f, nuse

def clim_series(Y0, Y1, minyears, stations_mask=None):
    stack = []
    for Y in range(Y0, Y1 + 1):
        a, b = window(Y)
        dates = np.arange(np.datetime64(f'{Y-1}-07-28'), np.datetime64(f'{Y}-09-01'))
        mm = dates.astype('datetime64[M]')
        keep = ~(((mm.astype(int) % 12) == 1) & ((dates - mm.astype('datetime64[D]')).astype(int) == 28))
        stack.append(Mf[sidx, a:b][:, keep])
    X = np.stack(stack)
    cnt = (~np.isnan(X)).sum(axis=0)
    m = np.nanmean(X, axis=0)
    m[cnt < minyears] = np.nan
    if stations_mask is not None:
        m[~stations_mask] = np.nan
    return m

ref0 = np.datetime64('2001-07-28')
RC0 = int((np.datetime64('2001-11-01') - ref0).astype(int))
RC1 = int((np.datetime64('2002-03-31') - ref0).astype(int))
RJ1 = int((np.datetime64('2002-01-01') - ref0).astype(int))

def run_clim(m):
    g, nuse, _ = interp(m)
    o, f = onset(g, cold=(RC0, RC1))
    return np.where(f <= 1, o - RJ1 + 1, -9999), f, nuse

if __name__ == '__main__':
    years = [y for y in list(range(1957, 1965)) + list(range(1974, 2021))]
    out = {}
    for Y in years:
        doy, f, nuse = run_year(Y)
        out[str(Y)] = (doy, f)
        print(Y, 'stations', nuse, 'flags', np.bincount(f, minlength=4).tolist(), flush=True)
    # eras: use only stations that have a climatology in both eras so the comparison is like-for-like
    e1 = clim_series(1957, 1964, 6); e2 = clim_series(2011, 2020, 8)
    both = (~np.isnan(e1)).mean(1) >= 0.95
    both &= (~np.isnan(e2)).mean(1) >= 0.95
    print('stations in both eras', both.sum())
    for key, m in [('e1', clim_series(1957, 1964, 6, both)), ('e2', clim_series(2011, 2020, 8, both)),
                   ('n9120', clim_series(1991, 2020, 20)), ('n8110', clim_series(1981, 2010, 20))]:
        doy, f, nuse = run_clim(m)
        out[key] = (doy, f)
        print(key, 'stations', nuse, 'flags', np.bincount(f, minlength=4).tolist(), flush=True)
    np.savez_compressed(os.path.join(BASE, 'grid_onsets.npz'),
                        keys=np.array(list(out.keys())),
                        doy=np.array([v[0] for v in out.values()]), flag=np.array([v[1] for v in out.values()]))
    json.dump({'era_stations': [ids[i] for i in sidx[both]]}, open(os.path.join(BASE, 'era_stations.json'), 'w'))
