#!/bin/bash
# Download every input into tools/work/ (GitHub only; about 3 GB of GSOD archives pass through, ~350 MB stay).
set -e
cd "$(dirname "$0")"
mkdir -p work/meta
# 1. CMA station list, 1 km DEM and boundary lines bundled with nmcdev/meteva (Apache-2.0)
if [ ! -d work/nmcdev_meteva ]; then
  git clone -q --depth 1 --filter=blob:none --sparse https://github.com/nmcdev/meteva work/nmcdev_meteva
  git -C work/nmcdev_meteva sparse-checkout set meteva/resources/stations meteva/resources/maps
fi
# 2. NOAA ISD station history, as packaged by rOpenSci GSODR (MIT)
[ -f work/isd_history.rda ] || curl -sfL -o work/isd_history.rda https://raw.githubusercontent.com/ropensci/GSODR/main/inst/extdata/isd_history.rda
# 3. Meteostat station metadata (CC BY 4.0) for stations missing from the two lists above
for s in 45004 45005 45007 45011 45032 45035 45039 45044 46686 46689 46692 46694 46695 46696 46697 46699 46706 46708 46730 46734 46735 46736 46740 46741 46743 46744 46745 46746 46748 46749 46750 46751 46752 46753 46754 46755 46756 46757 46758 46759 46760 46761 46762 46763 46764 46765 46766 46769 46770 46772 46777 46780 46810 51848 51895 52378 52441 53231 53546 53887 54587 54656 54868 54952 56075 56265 56294 56741 56945 57073 57451 57504 57515 57602 57622 58144 58239 58445 58473 58653 58726 58765 58849 58853 58965 58968 58969 58974 59158 59345 59358 59362 59559 59562 59567 59792 59985 59995 59997; do
  [ -f work/meta/$s.json ] || curl -sf -o work/meta/$s.json https://raw.githubusercontent.com/meteostat/weather-stations/master/stations/$s.json || true
done
# 4. NOAA GSOD daily data 1951-2020 from the CalebBell/gsod mirror; only China-block stations are kept
python3 fetch_gsod.py 1951 2020
