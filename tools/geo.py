"""Projection, China mask, analysis grid, DEM sampling and SVG path strings."""
import os, json, numpy as np, shapefile, netCDF4
from shapely.geometry import LineString, MultiLineString, Polygon, box, Point
from shapely.ops import polygonize, unary_union
from shapely import prepared, vectorized

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work')
MAPS = os.path.join(BASE, 'nmcdev_meteva/meteva/resources/maps/')
R = 6371.0
LON0, LAT1, LAT2, LAT0 = 105.0, 25.0, 47.0, 0.0
p1, p2, p0 = np.radians(LAT1), np.radians(LAT2), np.radians(LAT0)
n_ = (np.sin(p1) + np.sin(p2)) / 2
C_ = np.cos(p1) ** 2 + 2 * n_ * np.sin(p1)
rho0 = R * np.sqrt(C_ - 2 * n_ * np.sin(p0)) / n_

def proj(lon, lat):
    lon = np.asarray(lon, float); lat = np.asarray(lat, float)
    rho = R * np.sqrt(C_ - 2 * n_ * np.sin(np.radians(lat))) / n_
    th = n_ * np.radians(lon - LON0)
    return rho * np.sin(th), rho0 - rho * np.cos(th)  # km, y up

def inv(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    rho = np.sqrt(x ** 2 + (rho0 - y) ** 2)
    th = np.arctan2(x, rho0 - y)
    lat = np.degrees(np.arcsin((C_ - (rho * n_ / R) ** 2) / (2 * n_)))
    lon = LON0 + np.degrees(th / n_)
    return lon, lat

def read_lines(name):
    r = shapefile.Reader(os.path.join(MAPS, name))
    out = []
    for s in r.shapes():
        if not s.points:
            continue
        pts = s.points; parts = list(s.parts) + [len(pts)]
        for i in range(len(parts) - 1):
            seg = pts[parts[i]:parts[i + 1]]
            if len(seg) >= 2:
                out.append(np.array(seg))
    return out

if __name__ == '__main__':
    prov = read_lines('Province')
    nat = read_lines('NationalBorder')
    world = read_lines('worldmap')

    # China land polygon (lon/lat) from province rings
    polys = list(polygonize(unary_union([LineString(p) for p in prov])))  # node the rings first
    china_ll = unary_union([p.buffer(0) for p in polys])
    print('china polygons', len(polys), 'area deg2', round(china_ll.area, 1))

    def P(arr):
        x, y = proj(arr[:, 0], arr[:, 1]); return np.c_[x, y]

    def boxP(lon0, lat0, lon1, lat1):
        # densify the lon/lat box so its edges follow parallels/meridians after projection
        b = box(lon0, lat0, lon1, lat1).segmentize(0.25)
        return Polygon(P(np.array(b.exterior.coords)))

    china_xy = unary_union([Polygon(P(np.array(p.exterior.coords))).buffer(0) for p in polys])
    # main frame: exclude SCS south of 17N
    minx, miny, maxx, maxy = china_xy.intersection(boxP(70, 17.0, 140, 56)).bounds
    print('china xy bounds', minx, miny, maxx, maxy)

    # analysis grid (km)
    CELL = 20.0
    gx0, gy0 = np.floor(minx / CELL) * CELL - CELL, np.floor(miny / CELL) * CELL - CELL
    gx1, gy1 = np.ceil(maxx / CELL) * CELL + CELL, np.ceil(maxy / CELL) * CELL + CELL
    xs = np.arange(gx0 + CELL / 2, gx1, CELL); ys = np.arange(gy1 - CELL / 2, gy0, -CELL)  # rows top->bottom
    GX, GY = np.meshgrid(xs, ys)
    inside = vectorized.contains(china_xy, GX, GY)
    near = vectorized.contains(china_xy.buffer(CELL * 1.6), GX, GY)
    glon, glat = inv(GX, GY)
    print('grid', GX.shape, 'inside', inside.sum(), 'near', near.sum())

    # DEM block-mean at 0.1 deg then sample nearest
    d = netCDF4.Dataset(os.path.join(BASE, 'nmcdev_meteva/meteva/resources/stations/dem_0.00833.nc'))
    lat = d['lat'][:]; lon = d['lon'][:]
    Z = d['data0'][0, 0, 0, 0, :, :].astype(np.float32)
    Z = np.ma.filled(Z, np.nan) if np.ma.isMaskedArray(Z) else Z
    Z[Z < -500] = np.nan
    f = 12  # ~0.1 deg
    H, W = (Z.shape[0] // f) * f, (Z.shape[1] // f) * f
    Zb = np.nanmean(Z[:H, :W].reshape(H // f, f, W // f, f), axis=(1, 3))
    latb = lat[:H].reshape(-1, f).mean(1); lonb = lon[:W].reshape(-1, f).mean(1)
    print('dem', Z.shape, 'lat', float(lat[0]), float(lat[-1]), 'lon', float(lon[0]), float(lon[-1]))
    # cell mean over ~0.2deg: average 3x3 of the 0.1 deg blocks around the nearest block
    li = np.clip(np.round((glat - latb[0]) / (latb[1] - latb[0])).astype(int), 1, len(latb) - 2)
    lj = np.clip(np.round((glon - lonb[0]) / (lonb[1] - lonb[0])).astype(int), 1, len(lonb) - 2)
    acc = np.zeros_like(GX); cnt = np.zeros_like(GX)
    for di_ in (-1, 0, 1):
        for dj in (-1, 0, 1):
            v = Zb[li + di_, lj + dj]
            ok = ~np.isnan(v)
            acc[ok] += v[ok]; cnt[ok] += 1
    gz = np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0)
    gz = np.maximum(gz, 0)
    print('grid elev range', float(gz[inside].min()), float(gz[inside].max()))

    np.savez_compressed(os.path.join(BASE, 'grid.npz'), GX=GX, GY=GY, glon=glon, glat=glat, gz=gz,
                        inside=inside, near=near, cell=CELL, x0=xs[0], y0=ys[0])

    # ---- SVG geometry (pixel space) ----
    PADL, PADR, PADT, PADB = 40, 40, 30, 30
    WPX = 1000.0
    s = (WPX - PADL - PADR) / (maxx - minx)
    HPX = (maxy - miny) * s + PADT + PADB
    def to_px(x, y):
        return (x - minx) * s + PADL, (maxy - y) * s + PADT

    def path_from_lines(lines, simplify_km=2.5, clip=None, transform=None):
        out = []
        for arr in lines:
            ls = LineString(P(arr))
            if clip is not None:
                ls = ls.intersection(clip)
                if ls.is_empty:
                    continue
            ls = ls.simplify(simplify_km)
            geoms = getattr(ls, 'geoms', [ls])
            for g in geoms:
                if g.is_empty or g.geom_type not in ('LineString', 'LinearRing'):
                    continue
                c = np.array(g.coords)
                if transform:
                    px, py = transform(c[:, 0], c[:, 1])
                else:
                    px, py = to_px(c[:, 0], c[:, 1])
                out.append('M' + 'L'.join(f'{a:.1f},{b:.1f}' for a, b in zip(px, py)))
        return ''.join(out)

    frame_xy = box(minx - PADL / s, miny - PADB / s, maxx + PADR / s, maxy + PADT / s)

    def poly_path(geom, transform=None, simplify_km=2.0):
        geom = geom.simplify(simplify_km)
        out = []
        for pg in getattr(geom, 'geoms', [geom]):
            if pg.area < 30:
                continue
            for ring in [pg.exterior] + list(pg.interiors):
                c = np.array(ring.coords)
                px, py = (transform or to_px)(c[:, 0], c[:, 1])
                out.append('M' + 'L'.join(f'{a:.1f},{b:.1f}' for a, b in zip(px, py)) + 'Z')
        return ''.join(out)

    # world coastline/borders within frame (lines), skip parts that sit on China's border
    world_xy = []
    for arr in world:
        if arr[:, 0].max() < 55 or arr[:, 0].min() > 155 or arr[:, 1].max() < -5 or arr[:, 1].min() > 62:
            continue
        world_xy.append(arr)
    geo = {
        'W': WPX, 'H': round(HPX, 1),
        'china': poly_path(china_xy.intersection(frame_xy)),
        'prov': path_from_lines(prov, 3.0, clip=frame_xy),
        'nat': path_from_lines(nat, 2.0, clip=frame_xy.intersection(boxP(70, 17.0, 140, 56))),
        'world': path_from_lines(world_xy, 4.0, clip=frame_xy),
        'grid': {'x0': float(to_px(xs[0], ys[0])[0]), 'y0': float(to_px(xs[0], ys[0])[1]),
                 'dx': float(CELL * s), 'nx': int(GX.shape[1]), 'ny': int(GX.shape[0])},
        'scale_km_px': float(s),
    }
    # South China Sea inset (lower right): lon 106-123, lat 3-23.5
    ib_xy = boxP(106.5, 2.8, 123.0, 23.6)
    bx0, by0, bx1, by1 = ib_xy.bounds
    iw = 104.0
    si = iw / (bx1 - bx0)
    ih = (by1 - by0) * si
    ox, oy = WPX - iw - 6, HPX - ih - 6
    def to_inset(x, y):
        return (x - bx0) * si + ox, (by1 - y) * si + oy
    geo['inset'] = {
        'x': ox, 'y': oy, 'w': iw, 'h': round(ih, 1),
        'nat': path_from_lines(nat, 3.0, clip=ib_xy, transform=to_inset),
        'world': path_from_lines(world_xy, 5.0, clip=ib_xy, transform=to_inset),
        'china': poly_path(china_xy.intersection(ib_xy), transform=to_inset, simplify_km=3.0),
    }
    json.dump(geo, open(os.path.join(BASE, 'geo.json'), 'w'))
    json.dump({'minx': minx, 'maxy': maxy, 's': s, 'PADL': PADL, 'PADT': PADT}, open(os.path.join(BASE, 'pxform.json'), 'w'))
    print('svg H', HPX, 'sizes', {k: len(v) for k, v in geo.items() if isinstance(v, str)})
