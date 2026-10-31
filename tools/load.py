"""Load parsed GSOD China-block TSVs into a station x day matrix of daily mean temperature (degC)."""
import os, glob, numpy as np, pandas as pd

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work')
D0 = np.datetime64('1951-01-01')
D1 = np.datetime64('2020-12-31')
NDAY = int((D1 - D0).astype(int)) + 1

def load(min_cnt=4):
    frames = []
    for f in sorted(glob.glob(os.path.join(BASE, 'parsed', '*.tsv'))):
        if os.path.getsize(f) == 0:
            continue
        df = pd.read_csv(f, sep='\t', header=None,
                         names=['stn', 'wban', 'date', 'temp', 'cnt', 'tmax', 'tmin'],
                         dtype={'stn': str, 'wban': str, 'date': str})
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df = df[(df.temp < 9000) & df.stn.str.endswith('0')]
    df['cnt'] = pd.to_numeric(df.cnt, errors='coerce')
    stats = {'rows_all': len(df), 'rows_lowcnt': int((df.cnt < min_cnt).sum())}
    df = df[df.cnt >= min_cnt]
    df['t'] = (df.temp - 32.0) * 5.0 / 9.0
    d = pd.to_datetime(df.date, format='%Y%m%d').values.astype('datetime64[D]')
    df['di'] = (d - D0).astype(int)
    df['wmo'] = df.stn.str[:5]
    # a station may appear under >1 USAF/WBAN combo; average duplicates per wmo-day
    g = df.groupby(['wmo', 'di']).t.mean().reset_index()
    ids = sorted(g.wmo.unique())
    idx = {s: i for i, s in enumerate(ids)}
    M = np.full((len(ids), NDAY), np.nan, dtype=np.float32)
    M[g.wmo.map(idx).values, g.di.values] = g.t.values
    return ids, M, stats

if __name__ == '__main__':
    ids, M, stats = load()
    print(stats, len(ids), M.shape)
    np.savez_compressed(os.path.join(BASE, 'tmat.npz'), ids=np.array(ids), M=M)
