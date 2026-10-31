"""Station metadata: coordinates (CMA list > ISD history > Meteostat) and Chinese names."""
import os, json, glob, numpy as np, pandas as pd
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work')
ST = os.path.join(BASE, 'nmcdev_meteva/meteva/resources/stations/')

alt = {}
for l in open(ST + 'sta2411_alt.txt'):
    t = l.split()
    if len(t) >= 4 and t[0].isdigit() and len(t[0]) == 5:
        alt[t[0]] = (float(t[1]), float(t[2]), float(t[3]))
names = {}
for l in open(ST + 'station_id_pro_county.txt', encoding='utf-8'):
    t = l.split()
    if len(t) >= 3 and t[0].isdigit() and len(t[0]) == 5 and t[0] not in names:
        names[t[0]] = (t[1], t[-1])

isd = pd.read_csv(os.path.join(BASE, 'isd_cn.csv'), dtype=str)
isd['wmo'] = isd.STNID.str[:5]
isd = isd[isd.STNID.str[5] == '0'].copy()
isd['END'] = isd.END.astype(int)
isd = isd.sort_values('END', ascending=False)
isdm = {}
for _, r in isd.iterrows():
    if r.wmo in isdm or pd.isna(r.LAT):
        continue
    try:
        isdm[r.wmo] = (float(r.LON), float(r.LAT), float(r['ELEV(M)']), r.NAME)
    except ValueError:
        pass
ms = {}
for f in glob.glob(os.path.join(BASE, 'meta', '*.json')):
    j = json.load(open(f))
    loc = j['location']
    ms[j['id']] = (loc['longitude'], loc['latitude'], loc['elevation'] or 0, j['name'].get('en'))

# Chinese names for stations not in the CMA county list (HK/Macau/Taiwan/historical mainland numbers)
EXTRA = {
    '45004': ('香港', '京士柏'), '45005': ('香港', '香港天文台'), '45007': ('香港', '赤鱲角'), '45011': ('澳门', '澳门'),
    '46686': ('台湾', '桃园机场'), '46689': ('台湾', '马祖'), '46692': ('台湾', '台北'), '46694': ('台湾', '基隆'),
    '46695': ('台湾', '彭佳屿'), '46696': ('台湾', '台北松山'), '46697': ('台湾', '桃园'), '46699': ('台湾', '花莲'),
    '46706': ('台湾', '苏澳'), '46708': ('台湾', '宜兰'), '46730': ('台湾', '东吉岛'), '46734': ('台湾', '澎湖马公'),
    '46735': ('台湾', '澎湖'), '46736': ('福建', '金门'), '46740': ('台湾', '高雄机场'), '46741': ('台湾', '台南'),
    '46743': ('台湾', '台南机场'), '46744': ('台湾', '高雄'), '46745': ('台湾', '冈山'), '46746': ('台湾', '嘉义机场'),
    '46748': ('台湾', '嘉义'), '46749': ('台湾', '台中'), '46750': ('台湾', '屏东南'), '46751': ('台湾', '台中水湳'),
    '46752': ('台湾', '恒春'), '46753': ('台湾', '阿里山'), '46754': ('台湾', '大武'), '46755': ('台湾', '玉山'),
    '46756': ('台湾', '新竹机场'), '46757': ('台湾', '新竹'), '46758': ('台湾', '屏东北'), '46759': ('台湾', '恒春'),
    '46760': ('台湾', '台东志航'), '46761': ('台湾', '成功'), '46762': ('台湾', '兰屿'), '46763': ('台湾', '花莲机场'),
    '46764': ('台湾', '宜兰'), '46765': ('台湾', '日月潭'), '46766': ('台湾', '台东'), '46769': ('台湾', '龙潭'),
    '46770': ('台湾', '清泉岗'), '46772': ('台湾', '佳冬'), '46777': ('台湾', '梧栖'), '46810': ('广东', '东沙岛'),
    '58849': ('福建', '马祖'), '58965': ('台湾', '桃园'), '58968': ('台湾', '台北'), '58969': ('台湾', '台北松山'),
    '59158': ('台湾', '台中'), '59358': ('台湾', '台南'), '59362': ('台湾', '花莲'), '59559': ('台湾', '恒春'),
    '59562': ('台湾', '台东'),
    '51848': ('新疆', '安迪尔'), '52441': ('甘肃', '北山'), '53546': ('内蒙古', '新街'), '53887': ('山西', '长治'),
    '56075': ('甘肃', '郎木寺'), '56265': ('四川', '乾宁'), '56294': ('四川', '成都'), '56741': ('云南', '六库'),
    '56945': ('云南', '孟定'), '57073': ('河南', '洛阳'), '57504': ('四川', '内江'), '57515': ('重庆', '重庆沙坪坝'),
    '57602': ('四川', '泸州'), '58144': ('江苏', '清江'), '58445': ('浙江', '天目山'), '58726': ('福建', '七仙山'),
}
PROV_FIX = {'内蒙': '内蒙古', '黑龙': '黑龙江'}

def get(wmo):
    if wmo in alt:
        lon, lat, z = alt[wmo]
    elif wmo in isdm:
        lon, lat, z, _ = isdm[wmo]
    elif wmo in ms:
        lon, lat, z, _ = ms[wmo]
    else:
        return None
    if z < -100 or z > 9000:
        z = isdm.get(wmo, (0, 0, 0))[2]
    if wmo in names:
        prov, nm = names[wmo]
    elif wmo in EXTRA:
        prov, nm = EXTRA[wmo]
    else:
        prov, nm = None, None
    return dict(id=wmo, lon=lon, lat=lat, z=z, prov=PROV_FIX.get(prov, prov), name=nm)

if __name__ == '__main__':
    ids = list(np.load(os.path.join(BASE, 'tmat.npz'))['ids'])
    meta = {}
    for s in ids:
        m = get(s)
        if m and m['name'] and np.isfinite(m['z']) and np.isfinite(m['lat']) and np.isfinite(m['lon']):
            meta[s] = m
    print('with coords', len(meta), 'of', len(ids), 'unnamed', sum(1 for m in meta.values() if not m['name']))
    json.dump(meta, open(os.path.join(BASE, 'meta.json'), 'w'), ensure_ascii=False)
