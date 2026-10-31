"""Extract the China-block part (WMO 45, 46, 50-59) of NOAA's ISD station history from GSODR's isd_history.rda."""
import os, pyreadr
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work')
df = pyreadr.read_r(os.path.join(BASE, 'isd_history.rda'))['isd_history']
cn = df[df.STNID.str[:2].isin(['45', '46'] + [str(i) for i in range(50, 60)])]
cn.to_csv(os.path.join(BASE, 'isd_cn.csv'), index=False)
print(len(cn), 'stations')
