# 数据处理流程

网页本身不需要构建，`assets/data.js` 里已经是处理好的数据。只有想重新计算（比如补上 2021 年以后的数据、改插值方法）时才需要这里的脚本。

需要 Python 3.10+，以及：

```bash
pip install numpy pandas scipy shapely pyshp netCDF4 py7zr pyreadr ephem
```

## 输入（由 `fetch_inputs.sh` 下载到 `tools/work/`）

| 输入 | 来源 |
|---|---|
| 逐日平均气温 | NOAA 全球地面日值 GSOD，取自 GitHub 镜像 [CalebBell/gsod](https://github.com/CalebBell/gsod)（每年一个 7z，1929–2020）。只保留 WMO 区号 45、46、50–59 的站点（中国大陆、港澳台） |
| 国家站中文名、坐标、海拔；1 公里地形；国界、省界线划 | [nmcdev/meteva](https://github.com/nmcdev/meteva) 的 `meteva/resources/` |
| 港澳台及已撤销老站号的位置 | NOAA ISD 站点历史（经 [ropensci/GSODR](https://github.com/ropensci/GSODR) 打包的 `isd_history.rda`）和 [Meteostat 站点表](https://github.com/meteostat/weather-stations) |

## 运行顺序（在 `tools/` 目录下）

```bash
./fetch_inputs.sh        # 下载上面的输入，并从 GSOD 年度包中抽出中国站点 → work/parsed/*.tsv
python3 isd.py           # ISD 站点历史 → work/isd_cn.csv
python3 load.py          # 站点 × 日期的日平均气温矩阵 → work/tmat.npz
python3 stations.py      # 每站每年的入春日，以及各时期“常年”入春日 → work/station_onsets.json
python3 meta.py          # 站点中文名、坐标、海拔 → work/meta.json
python3 geo.py           # 投影、20 公里网格、地形采样、SVG 底图 → work/grid.npz、work/geo.json
python3 grid_onset.py    # 逐日气温插值到网格，再逐格判定入春 → work/grid_onsets.npz（约 2 分钟）
python3 jieqi.py         # 1956–2026 年二十四节气日期（北京时间） → work/jieqi.json
python3 build_data.py    # 打包成 ../assets/data.js
```

## 各文件的作用

- `onset.py`：入春判定，依据 GB/T 42074—2022《气候季节划分》。滑动平均取当天及前 4 天；滑动平均连续 5 天 ≥10℃ 时，在对应 9 天的实测日平均气温里取第一个 ≥10℃ 的日期。入春、入冬都做“二次判断”（之后若出现连续 5 天滑动平均不达标，比较两段之间达标与不达标的天数）。冬半年没有冬季过程记为无冬，到 8 月底仍未入春记为无春。
- `stations.py`：缺测 3 天以内线性插补；某段资料缺失会影响判定时记为无资料。“常年”按时期的平均气温曲线计算（1957–1964、2011–2020、1981–2010、1991–2020）。
- `grid_onset.py`：每天先对全部站点做“气温 ~ 纬度 + 海拔”回归，残差按最近 10 站反距离平方插值，再对每个格点判定入春。两代对比只用两个时期资料都完整的站点。
- `geo.py`：Albers 等积圆锥投影（中央经线 105°E，标准纬线 25°N、47°N），南海诸岛放在右下角插图。
- `build_data.py`：网格压缩成 gzip + base64，连同站点、节气、底图路径写进 `window.SPRING`。

## 补 2021 年以后的数据

把 NOAA 官网的 GSOD 年度包（https://www.ncei.noaa.gov/data/global-summary-of-the-day/archive/ 下的 `YYYY.tar.gz`）解开，按 `fetch_gsod.py` 的格式写成 `work/parsed/YYYY.tsv`（列：站号、WBAN、日期、日平均气温 ℉、观测次数、最高、最低），再把各脚本里的 `2020` 改成新的年份，从 `load.py` 开始重跑。
