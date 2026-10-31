"""Pack stations, grids, solar terms and map geometry into one JS data file for the page."""
import os, sys, json, gzip, base64, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from geo import proj
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'assets')
os.makedirs(OUT, exist_ok=True)

geo = json.load(open(os.path.join(BASE, 'geo.json')))
px = json.load(open(os.path.join(BASE, 'pxform.json')))
meta = json.load(open(os.path.join(BASE, 'meta.json')))
so = json.load(open(os.path.join(BASE, 'station_onsets.json')))
jq = json.load(open(os.path.join(BASE, 'jieqi.json')))
G = np.load(os.path.join(BASE, 'grid.npz'))
R = np.load(os.path.join(BASE, 'grid_onsets.npz'))
era_st = set(json.load(open(os.path.join(BASE, 'era_stations.json')))['era_stations'])

def topx(lon, lat):
    x, y = proj(lon, lat)
    return (x - px['minx']) * px['s'] + px['PADL'], (px['maxy'] - y) * px['s'] + px['PADT']

# ---------------- stations ----------------
years = so['years']
y0 = years.index(1957)
YRS = years[y0:]
doy = np.array(so['doy'])[:, y0:]
fl = np.array(so['flag'])[:, y0:]
stations = []
for i, sid in enumerate(so['ids']):
    if sid not in meta:
        continue
    good = int((fl[i] != 3).sum())
    if good < 8:
        continue
    m = meta[sid]
    x, y = topx(m['lon'], m['lat'])
    o = [int(v) if f != 3 else None for v, f in zip(doy[i], fl[i])]
    c = {k: [int(v[0][i]), int(v[1][i])] for k, v in so['clim'].items()}
    stations.append({'id': sid, 'n': m['name'], 'p': m['prov'], 'z': round(m['z']), 'lat': round(m['lat'], 2),
                     'lon': round(m['lon'], 2), 'x': round(x, 1), 'y': round(y, 1),
                     'o': o, 'f': ''.join(str(v) for v in fl[i]), 'c': c, 'era': sid in era_st})
print('stations', len(stations))

# headline: station-level era difference (stations with normal onset in both era climatologies)
d = []
for s in stations:
    e1, e2 = s['c']['e1'], s['c']['e2']
    if e1[1] == 0 and e2[1] == 0:
        d.append(e2[0] - e1[0])
d = np.array(d)
print('era diff stations', len(d), 'median', np.median(d), 'mean', d.mean().round(1), 'earlier share', (d < 0).mean().round(2),
      'p10/p90', np.percentile(d, [10, 90]))

# ---------------- grids ----------------
near = G['near']; inside = G['inside']
NC = int(near.sum())
inside_near = inside[near]
keys = [str(k) for k in R['keys']]
blobs = []
offsets = {}
pos = 0
def bits(b):
    return np.packbits(b.astype(np.uint8), bitorder='little').tobytes()
mask_blob = bits(near.ravel()) + bits(inside.ravel())
blobs.append(mask_blob); pos += len(mask_blob)
for k, dd, ff in zip(keys, R['doy'], R['flag']):
    code = np.zeros(NC, np.uint8)
    normal = ff <= 1
    code[normal] = np.clip(dd[normal], -27, 225) + 28
    code[ff == 2] = 254
    code[ff == 3] = 0
    b = code.tobytes() + bits(ff == 1)
    offsets[k] = pos
    blobs.append(b); pos += len(b)
raw = b''.join(blobs)
gz = gzip.compress(raw, 9)
print('grid raw', len(raw), 'gz', len(gz))

# area curves for the headline: share of land in spring by day for e1/e2
# (computed in the page too; here only to print a sanity number)
def share(k, day):
    i = keys.index(k)
    dd, ff = R['doy'][i][inside_near], R['flag'][i][inside_near]
    ok = ff <= 1
    return float(((dd <= day) & ok).sum() / len(dd))
for day in (64, 79, 95, 110, 126):
    print('day', day, 'e1', round(share('e1', day), 3), 'e2', round(share('e2', day), 3))

# ---------------- solar terms ----------------
import datetime as dt
TERMS = ['冬至', '小寒', '大寒', '立春', '雨水', '惊蛰', '春分', '清明', '谷雨', '立夏', '小满', '芒种', '夏至']
jqd = {}
for Y, row in jq.items():
    Y = int(Y)
    jan1 = dt.date(Y, 1, 1)
    jqd[Y] = [(dt.date.fromisoformat(row[t]) - jan1).days + 1 for t in TERMS]
mean_jq = np.round(np.mean([jqd[Y] for Y in range(1957, 2021)], axis=0)).astype(int).tolist()

data = {
    'geo': geo, 'nx': int(near.shape[1]), 'ny': int(near.shape[0]), 'NC': NC,
    'gridKeys': keys, 'gridOff': offsets, 'grid': base64.b64encode(gz).decode(),
    'years': YRS, 'stations': stations,
    'terms': TERMS, 'jq': {str(k): v for k, v in jqd.items()}, 'jqMean': mean_jq,
    'eraN': len(era_st),
}
js = 'window.SPRING=' + json.dumps(data, ensure_ascii=False, separators=(',', ':')) + ';'
open(os.path.join(OUT, 'data.js'), 'w').write(js)
print('data.js', round(len(js) / 1e6, 2), 'MB', 'mean jieqi', dict(zip(TERMS, mean_jq)))
