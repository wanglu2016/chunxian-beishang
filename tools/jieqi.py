import os, ephem, json, math, datetime as dt
TERMS = [('冬至',270),('小寒',285),('大寒',300),('立春',315),('雨水',330),('惊蛰',345),('春分',0),('清明',15),
         ('谷雨',30),('立夏',45),('小满',60),('芒种',75),('夏至',90)]
def sunlon(d):
    s = ephem.Sun(); s.compute(d, epoch=d)
    return math.degrees(float(ephem.Ecliptic(s, epoch=d).lon)) % 360
def find(target, start):
    # step daily until crossing, then bisect
    d = ephem.Date(start)
    prev = (sunlon(d) - target + 180) % 360 - 180
    for _ in range(400):
        d2 = ephem.Date(d + 1)
        cur = (sunlon(d2) - target + 180) % 360 - 180
        if prev < 0 <= cur:
            lo, hi = d, d2
            for _ in range(40):
                mid = ephem.Date((lo + hi) / 2)
                if ((sunlon(mid) - target + 180) % 360 - 180) < 0: lo = mid
                else: hi = mid
            return ephem.Date(hi)
        d, prev = d2, cur
    raise RuntimeError
out = {}
for Y in range(1956, 2027):
    row = {}
    for name, lon in TERMS:
        start = f'{Y-1}/12/10' if name == '冬至' else f'{Y}/01/01'
        t = find(lon, start)
        bj = ephem.Date(t + 8 * ephem.hour).datetime().date()
        row[name] = bj.isoformat()
    out[Y] = row
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work', 'jieqi.json'), 'w'), ensure_ascii=False)
print(out[2024]); print(out[1960])
