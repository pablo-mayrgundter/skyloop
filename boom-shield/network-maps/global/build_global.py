#!/usr/bin/env python3
"""
BoomShield global corridor study (screening level).

Question: which long overland spines elsewhere in the world have economics like, or better than, the US
LAS-DEN-OMA-CHI-CLE-NYC spine (coverage.md / highmach.md)?

Pipeline (each stage cached in .cache/ so a container restart resumes):
  1. Restricted-airspace grid (Natural Earth 1:50m land + 40 km coastal buffer, routes.md convention), 0.25 deg.
  2. Demand per system: US = demand/.cache/inputs.pkl (DOT Table 1 + gravity fill, read-only); elsewhere a
     gravity model T_ij = k m_i m_j d^-beta (1 - exp(-(d/d0)^2)) x affinity, domestic and international
     legs separately, calibrated per country to national totals + published route anchors (evidence.csv,
     unverified search summaries; [U] = author estimates; [K] = background knowledge).
  3. Baselines per pair: subsonic M0.85; Boomless (same design aircraft, no shield): M1.15 over restricted
     ground, design Mach over open-water runs, including one-waypoint offshore tracks (coverage.md idea).
  4. Spines: hand-built corridors (the user's list) + a greedy builder per region (trunk.md style: seed with
     the densest pair's path on a Delaunay graph, grow by best marginal pax-h per overland km, stop below
     0.5x the average).  Routing = trunk.md: same aircraft, feeder M0.85/M1.15 to any ramp, t_trunk().
  5. Economics with highmach.md conventions: design A M2.2 (and M2.0 adaptive), cost base / integrated
     (aggressive x integration ratio), $86/h x 0.95, induced/day-return uplift from highmach's spine,
     fuel and CO2 netted with hm_aircraft fuel model.  Only restricted (overland) km pay for shield.
  6. Mach-cutoff day share (parametric climatology: latitude/elevation temperature, jet streams),
     jet-level wind exposure, terrain, regulation and politics flags.
  7. Outputs: CSVs + figures/global_map.png.

Usage: python3 build_global.py [--stage N]     one process; ~20-40 min cold.
"""
from __future__ import annotations

import argparse
import math
import pickle
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, cKDTree

HERE = Path(__file__).resolve().parent
NM = HERE.parent
BS = NM.parent
CACHE = HERE / ".cache"
FIG = HERE / "figures"
CACHE.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)
sys.dont_write_bytecode = True
sys.path.insert(0, str(NM / "trunk"))
sys.path.insert(0, str(NM / "highmach"))
import build_trunk as bt  # noqa: E402  (time model, read-only)
import hm_aircraft as ha  # noqa: E402  (fuel model, read-only)

R_E = 6371.0
RES = 0.25
CUT_MIN = 10.0
RANGE_KM = 7000.0      # nonstop cap on pair great circle: Overture-class 4,250 nmi (7,870 km) less detour [U];
                       # highmach concept A is sized for 5,000 km, run as a sensitivity (RANGE_SENS)
RANGE_SENS = 5000.0
VOT, AVAIL = 86.0, 0.95
FUEL_USD_KG, CO2_KG, SCC = 0.80, 3.16, 190.0
FUEL_SUB_PAX_KM = 0.090 / CO2_KG
LF = 0.75
SEATS = 65


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def cached(name, fn):
    f = CACHE / f"{name}.pkl"
    if f.exists():
        return pickle.loads(f.read_bytes())
    out = fn()
    f.write_bytes(pickle.dumps(out))
    return out


# =====================================================================================
# 0. Geometry
# =====================================================================================
def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def gc(lat1, lon1, lat2, lon2):
    a, b = xyz(lat1, lon1), xyz(lat2, lon2)
    return R_E * 2 * np.arcsin(np.clip(np.linalg.norm(a - b, axis=-1) / 2, 0, 1))


def slerp(p, q, f):
    """Points at fractions f along the great circle p->q. p, q = (lat, lon)."""
    a, b = xyz(*p), xyz(*q)
    om = math.acos(float(np.clip(a @ b, -1, 1)))
    f = np.asarray(f, float)
    if om < 1e-9:
        return np.full(f.shape, p[0]), np.full(f.shape, p[1])
    v = (np.sin((1 - f) * om)[:, None] * a + np.sin(f * om)[:, None] * b) / math.sin(om)
    return np.degrees(np.arcsin(np.clip(v[:, 2], -1, 1))), np.degrees(np.arctan2(v[:, 1], v[:, 0]))


def bearing(lat1, lon1, lat2, lon2):
    p1, p2, dl = np.radians(lat1), np.radians(lat2), np.radians(lon2 - lon1)
    x = np.sin(dl) * np.cos(p2)
    y = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return (np.degrees(np.arctan2(x, y)) + 360) % 360


def restricted_grid():
    """uint8 grid [lat, lon]: 1 = restricted (land or within 40 km of coast)."""
    f = CACHE / "restricted_grid_025.npy"
    if f.exists():
        return np.load(f)
    sys.path.insert(0, str(NM / "routes"))
    import build_routes as br
    lm = br.LandMask(40.0)
    lat = -90 + RES / 2 + RES * np.arange(int(180 / RES))
    lon = -180 + RES / 2 + RES * np.arange(int(360 / RES))
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    la, lo = LA.ravel(), LO.ravel()
    out = np.zeros(la.size, np.uint8)
    for s in range(0, la.size, 200000):
        _, r = lm.classify(la[s:s + 200000], lo[s:s + 200000])
        out[s:s + 200000] = r
        log(f"  grid {s / la.size:.0%}")
    g = out.reshape(LA.shape)
    np.save(f, g)
    return g


GRID = None


def restricted(lat, lon):
    i = np.clip(((np.asarray(lat) + 90) / RES).astype(int), 0, GRID.shape[0] - 1)
    j = ((np.asarray(lon) + 180) / RES).astype(int) % GRID.shape[1]
    return GRID[i, j].astype(bool)


# Coarse plateau elevation boxes, km [K]: (lat0, lat1, lon0, lon1, elev_km, name)
PLATEAUS = [
    (28, 36.5, 78, 101, 4.5, "Tibetan Plateau"), (36, 39.5, 89, 99, 3.0, "Qaidam/Qilian"),
    (36.5, 41, 70, 80, 3.5, "Pamir/Karakoram"), (40, 44, 75, 86, 2.0, "Tien Shan"),
    (26, 34, 100, 108, 1.5, "Yunnan-Guizhou"), (34, 39, 103, 112, 1.3, "Loess Plateau"),
    (41, 50, 95, 116, 1.3, "Mongolian Plateau"), (27, 38, 45, 62, 1.2, "Iranian Plateau"),
    (31, 37, 62, 70, 1.8, "Afghan highlands"), (37, 41, 30, 44, 1.0, "Anatolia"),
    (5, 14, 36, 41, 2.2, "Ethiopian Highlands"), (-5, 2, 33, 38, 1.5, "East African Plateau"),
    (-29, -23, 25, 31, 1.5, "Highveld"), (-18, -8, 22, 34, 1.1, "Zambian plateau"),
    (12, 22, 74, 80, 0.6, "Deccan"), (-22, -14, -70, -65, 3.8, "Altiplano"),
    (-35, -22, -71, -68, 3.0, "Southern Andes crest"), (-14, 2, -79, -74, 2.8, "Northern Andes"),
    (2, 8, -77, -72, 2.0, "Colombian Andes"), (-20, -10, -50, -43, 0.9, "Brazilian Highlands"),
    (35, 45, -112, -105, 2.0, "Rockies/Colorado Plateau"), (37, 42.5, -120, -112, 1.5, "Great Basin"),
    (35, 45, -105, -100, 1.2, "High Plains"), (18, 25, -104, -97, 2.0, "Mexican Plateau"),
    (49, 54, -123, -114, 1.5, "Canadian Rockies"), (45.5, 47.5, 6, 14, 1.2, "Alps"),
    (21, 30, 39, 46, 0.8, "Najd"),
]


def elevation(lat, lon):
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    e = np.full(lat.shape, 0.3)
    for a0, a1, o0, o1, z, _ in PLATEAUS:
        m = (lat >= a0) & (lat <= a1) & (lon >= o0) & (lon <= o1)
        e = np.where(m, np.maximum(e, z), e)
    return e


def plateau_names(lat, lon):
    out = set()
    for a0, a1, o0, o1, z, nm in PLATEAUS:
        if z >= 1.0 and np.any((lat >= a0) & (lat <= a1) & (lon >= o0) & (lon <= o1)):
            out.add(f"{nm} ~{z:.1f} km")
    return sorted(out)


# =====================================================================================
# 1. Metros.  code, name, lat, lon, country, metro airport pax 2024 (M) [U], transfer share [U],
#    international share of local O&D [U].  local O&D mass = pax x (1 - transfer).
# =====================================================================================
OW = """
LON London 51.47 -0.45 GB 181 .15 .9
PAR Paris 49.0 2.55 FR 103 .25 .85
FRA Frankfurt 50.03 8.56 DE 62 .5 .85
MUC Munich 48.35 11.79 DE 41 .35 .75
BER Berlin 52.37 13.5 DE 25 .05 .8
HAM Hamburg 53.63 9.99 DE 14 0 .65
DUS Dusseldorf/Cologne 51.29 6.77 DE 33 .05 .85
AMS Amsterdam 52.31 4.76 NL 67 .35 1
BRU Brussels 50.9 4.48 BE 22 .1 1
ZRH Zurich 47.46 8.55 CH 31 .3 .97
GVA Geneva 46.24 6.11 CH 18 0 .95
VIE Vienna 48.11 16.57 AT 32 .3 .98
MAD Madrid 40.47 -3.57 ES 66 .25 .7
BCN Barcelona 41.3 2.08 ES 55 .05 .75
LIS Lisbon 38.77 -9.13 PT 35 .2 .85
ROM Rome 41.8 12.25 IT 50 .15 .65
MIL Milan 45.63 8.72 IT 50 .05 .75
ATH Athens 37.94 23.94 GR 32 .1 .7
IST Istanbul 41.26 28.74 TR 121 .35 .6
ANK Ankara 40.13 32.99 TR 14 0 .15
AYT Antalya 36.9 30.8 TR 39 0 .85
WAW Warsaw 52.17 20.97 PL 21 .1 .85
PRG Prague 50.1 14.26 CZ 14 0 1
BUD Budapest 47.44 19.26 HU 17 0 1
OTP Bucharest 44.57 26.08 RO 16 0 .85
SOF Sofia 42.7 23.4 BG 8 0 .9
BEG Belgrade 44.82 20.29 RS 8 .1 1
CPH Copenhagen 55.62 12.66 DK 30 .2 .95
ARN Stockholm 59.65 17.92 SE 27 .1 .75
OSL Oslo 60.19 11.1 NO 29 .1 .55
HEL Helsinki 60.32 24.95 FI 16 .3 .85
DUB Dublin 53.42 -6.27 IE 33 .05 1
MAN Manchester 53.35 -2.27 GB 31 0 .9
EDI Edinburgh 55.95 -3.37 GB 15 0 .6
MOW Moscow 55.97 37.41 RU 100 .25 .2
LED St_Petersburg 59.8 30.26 RU 22 0 .12
KZN Kazan 55.6 49.28 RU 5 0 .2
KUF Samara 53.5 50.16 RU 4 0 .15
UFA Ufa 54.56 55.87 RU 4 0 .2
SVX Yekaterinburg 56.74 60.8 RU 8 0 .25
TJM Tyumen 57.19 65.32 RU 2 0 .05
SGC Surgut 61.34 73.4 RU 2 0 .05
OMS Omsk 54.97 73.31 RU 1.6 0 .1
OVB Novosibirsk 55.01 82.65 RU 8 0 .2
KJA Krasnoyarsk 56.17 92.49 RU 4 0 .1
IKT Irkutsk 52.27 104.39 RU 3 0 .15
YKS Yakutsk 62.09 129.77 RU 1.5 0 0
KHV Khabarovsk 48.53 135.19 RU 3 0 .1
VVO Vladivostok 43.4 132.15 RU 3 0 .2
AER Sochi 43.45 39.96 RU 10 0 .05
MRV Mineralnye_Vody 44.22 43.08 RU 3 0 .05
KGD Kaliningrad 54.89 20.59 RU 3 0 0
GYD Baku 40.47 50.05 AZ 6 .1 1
TBS Tbilisi 41.67 44.95 GE 5 0 1
EVN Yerevan 40.15 44.4 AM 4 0 1
ALA Almaty 43.35 77.04 KZ 11 .05 .45
NQZ Astana 51.02 71.47 KZ 7 .05 .35
TAS Tashkent 41.26 69.28 UZ 9 .05 .85
FRU Bishkek 43.06 74.48 KG 3 0 .8
DXB Dubai 25.25 55.36 AE 92 .55 1
AUH Abu_Dhabi 24.43 54.65 AE 29 .45 1
SHJ Sharjah 25.33 55.52 AE 17 .05 1
DOH Doha 25.27 51.61 QA 52 .7 1
BAH Bahrain 26.27 50.63 BH 9 .15 1
KWI Kuwait 29.24 47.97 KW 15 .05 1
RUH Riyadh 24.96 46.7 SA 37 .05 .4
JED Jeddah 21.68 39.16 SA 49 .05 .5
DMM Dammam 26.47 49.8 SA 12 0 .45
MED Madinah 24.55 39.7 SA 10 0 .5
MCT Muscat 23.59 58.28 OM 15 .15 .9
AMM Amman 31.72 35.99 JO 9 .1 1
BEY Beirut 33.82 35.49 LB 7 0 1
TLV Tel_Aviv 32.01 34.89 IL 21 0 1
CAI Cairo 30.12 31.41 EG 30 .15 .85
BGW Baghdad 33.26 44.23 IQ 6 0 .8
THR Tehran 35.42 51.15 IR 25 0 .3
MHD Mashhad 36.24 59.64 IR 15 0 .15
SYZ Shiraz 29.54 52.59 IR 6 0 .15
DEL Delhi 28.56 77.1 IN 78 .1 .25
BOM Mumbai 19.09 72.87 IN 55 .05 .3
BLR Bengaluru 13.2 77.71 IN 41 .05 .15
MAA Chennai 12.99 80.17 IN 22 0 .2
HYD Hyderabad 17.24 78.43 IN 28 0 .15
CCU Kolkata 22.65 88.45 IN 21 0 .12
AMD Ahmedabad 23.07 72.63 IN 13 0 .2
COK Kochi 10.15 76.4 IN 11 0 .45
PNQ Pune 18.58 73.92 IN 10 0 .05
GOI Goa 15.38 73.83 IN 10 0 .1
JAI Jaipur 26.82 75.81 IN 6 0 .1
LKO Lucknow 26.76 80.89 IN 6 0 .2
GAU Guwahati 26.11 91.59 IN 6 0 .02
PAT Patna 25.59 85.09 IN 4 0 0
TRV Thiruvananthapuram 8.48 76.92 IN 5 0 .5
SXR Srinagar 33.99 74.77 IN 4 0 0
KHI Karachi 24.91 67.16 PK 9 0 .4
LHE Lahore 31.52 74.4 PK 6 0 .5
ISB Islamabad 33.55 72.83 PK 5 0 .6
KTM Kathmandu 27.7 85.36 NP 8 0 .6
DAC Dhaka 23.84 90.4 BD 10 0 .6
CMB Colombo 7.18 79.88 LK 10 .1 1
BKK Bangkok 13.69 100.75 TH 94 .12 .55
CNX Chiang_Mai 18.77 98.96 TH 10 0 .2
HKT Phuket 8.11 98.31 TH 16 0 .55
SIN Singapore 1.36 103.99 SG 68 .3 1
KUL Kuala_Lumpur 2.75 101.71 MY 57 .2 .5
PEN Penang 5.3 100.28 MY 8 0 .4
CGK Jakarta -6.13 106.66 ID 55 .05 .25
DPS Bali -8.75 115.17 ID 24 0 .5
SUB Surabaya -7.38 112.79 ID 20 0 .1
MNL Manila 14.51 121.02 PH 50 .05 .4
SGN Ho_Chi_Minh 10.82 106.66 VN 41 0 .3
HAN Hanoi 21.22 105.81 VN 29 0 .3
DAD Da_Nang 16.04 108.2 VN 15 0 .4
RGN Yangon 16.91 96.13 MM 6 0 .4
PNH Phnom_Penh 11.55 104.84 KH 5 0 .8
VTE Vientiane 17.99 102.56 LA 2.5 0 .6
PEK Beijing 40.08 116.58 CN 117 .15 .12
SHA Shanghai 31.19 121.6 CN 125 .1 .2
CAN Guangzhou 23.39 113.3 CN 76 .15 .15
SZX Shenzhen 22.64 113.81 CN 62 .05 .08
CTU Chengdu 30.58 103.95 CN 80 .05 .05
CKG Chongqing 29.72 106.64 CN 49 .05 .04
KMG Kunming 25.1 102.93 CN 47 .05 .06
XIY Xian 34.45 108.75 CN 47 .05 .04
HGH Hangzhou 30.23 120.43 CN 48 0 .07
WUH Wuhan 30.78 114.21 CN 32 0 .05
CSX Changsha 28.19 113.22 CN 30 0 .05
NKG Nanjing 31.74 118.86 CN 33 0 .05
XMN Xiamen 24.54 118.13 CN 28 0 .1
TAO Qingdao 36.36 120.09 CN 30 0 .06
SHE Shenyang 41.64 123.48 CN 23 0 .05
HRB Harbin 45.62 126.25 CN 23 0 .03
CGQ Changchun 43.99 125.68 CN 15 0 .02
DLC Dalian 38.97 121.54 CN 20 0 .05
URC Urumqi 43.9 87.47 CN 28 0 .03
KHG Kashgar 39.54 76.02 CN 4 0 0
LHW Lanzhou 36.52 103.62 CN 16 0 .01
XNN Xining 36.53 102.04 CN 10 0 0
INC Yinchuan 38.32 106.39 CN 12 0 .01
HET Hohhot 40.85 111.82 CN 13 0 .01
LXA Lhasa 29.3 90.91 CN 6 0 0
KWE Guiyang 26.54 106.8 CN 25 0 .02
NNG Nanning 22.61 108.17 CN 17 0 .03
HAK Haikou 19.93 110.46 CN 26 0 .04
SYX Sanya 18.3 109.41 CN 22 0 .04
TSN Tianjin 39.12 117.35 CN 20 0 .04
CGO Zhengzhou 34.52 113.84 CN 31 0 .04
TNA Jinan 36.86 117.22 CN 20 0 .03
TYN Taiyuan 37.75 112.63 CN 15 0 .02
FOC Fuzhou 25.94 119.66 CN 18 0 .06
HFE Hefei 31.99 116.97 CN 13 0 .02
NGB Ningbo 29.83 121.46 CN 14 0 .05
WNZ Wenzhou 27.91 120.85 CN 13 0 .03
HKG Hong_Kong 22.31 113.91 HK 53 .3 1
MFM Macau 22.15 113.59 MO 7 0 1
TPE Taipei 25.08 121.23 TW 45 .15 1
SEL Seoul 37.46 126.44 KR 96 .1 .75
TYO Tokyo 35.55 139.78 JP 125 .05 .35
OSA Osaka 34.43 135.24 JP 45 0 .5
FUK Fukuoka 33.59 130.45 JP 25 0 .3
CTS Sapporo 42.78 141.69 JP 24 0 .15
ULN Ulaanbaatar 47.65 106.82 MN 2 0 .8
JNB Johannesburg -26.14 28.25 ZA 21 .15 .35
CPT Cape_Town -33.97 18.6 ZA 11 0 .25
DUR Durban -29.61 31.12 ZA 6 0 .05
NBO Nairobi -1.32 36.93 KE 10 .3 .6
MBA Mombasa -4.03 39.6 KE 2 0 .2
ADD Addis_Ababa 8.98 38.8 ET 13 .6 .8
LOS Lagos 6.58 3.32 NG 10 0 .5
ABV Abuja 9.0 7.26 NG 4 0 .15
ACC Accra 5.61 -0.17 GH 3.5 0 .8
ABJ Abidjan 5.26 -3.93 CI 2.5 0 .9
DSS Dakar 14.67 -17.07 SN 3 0 1
DAR Dar_es_Salaam -6.88 39.2 TZ 3 0 .5
EBB Entebbe 0.04 32.44 UG 2 0 1
KGL Kigali -1.97 30.14 RW 1.2 .3 1
LAD Luanda -8.86 13.23 AO 4 0 .4
FIH Kinshasa -4.39 15.44 CD 1.5 0 .5
LUN Lusaka -15.33 28.45 ZM 2 0 .7
HRE Harare -17.93 31.09 ZW 1.3 0 .8
CMN Casablanca 33.37 -7.59 MA 11 0 .95
ALG Algiers 36.69 3.22 DZ 9 0 .6
TUN Tunis 36.85 10.23 TN 7 0 1
"""

CA = """
YYZ Toronto 43.68 -79.63 CA 50 .2 .55
YVR Vancouver 49.19 -123.18 CA 26 .1 .45
YUL Montreal 45.47 -73.74 CA 22 .1 .6
YYC Calgary 51.13 -114.01 CA 19 .05 .25
YEG Edmonton 53.31 -113.58 CA 8.4 0 .15
YOW Ottawa 45.32 -75.67 CA 5 0 .2
YWG Winnipeg 49.91 -97.24 CA 4.6 0 .15
YHZ Halifax 44.88 -63.51 CA 4.3 0 .1
YQB Quebec 46.79 -71.39 CA 1.8 0 .3
YXE Saskatoon 52.17 -106.7 CA 1.5 0 .1
YQR Regina 50.43 -104.67 CA 1.3 0 .1
YYJ Victoria 48.65 -123.43 CA 2.0 0 .05
YYT St_Johns 47.62 -52.75 CA 1.6 0 .05
YQT Thunder_Bay 48.37 -89.32 CA 0.8 0 0
"""

LA = """
SAO Sao_Paulo -23.43 -46.47 BR 70 .2 .2
RIO Rio_de_Janeiro -22.81 -43.25 BR 25 .05 .2
BSB Brasilia -15.87 -47.92 BR 15 .35 .02
CNF Belo_Horizonte -19.63 -43.97 BR 12 .05 .05
POA Porto_Alegre -29.99 -51.17 BR 7 0 .05
CWB Curitiba -25.53 -49.18 BR 7 0 .02
FLN Florianopolis -27.67 -48.55 BR 5 0 .1
SSA Salvador -12.91 -38.33 BR 8 0 .05
REC Recife -8.13 -34.92 BR 10 .1 .05
FOR Fortaleza -3.78 -38.53 BR 7 0 .05
NAT Natal -5.77 -35.37 BR 2.5 0 .05
BEL Belem -1.38 -48.48 BR 4 0 .02
MAO Manaus -3.04 -60.05 BR 3.5 0 .03
GYN Goiania -16.63 -49.22 BR 3 0 0
CGB Cuiaba -15.65 -56.12 BR 3 0 0
VIX Vitoria -20.26 -40.29 BR 3 0 0
BUE Buenos_Aires -34.82 -58.54 AR 25 .05 .45
COR Cordoba -31.32 -64.21 AR 4 0 .1
MDZ Mendoza -32.83 -68.79 AR 2.5 0 .15
SCL Santiago -33.39 -70.79 CL 26 .1 .4
LIM Lima -12.02 -77.11 PE 25 .15 .45
CUZ Cusco -13.54 -71.94 PE 5 0 0
BOG Bogota 4.7 -74.15 CO 45 .1 .3
MDE Medellin 6.16 -75.42 CO 13 0 .15
CLO Cali 3.54 -76.38 CO 6 0 .1
CTG Cartagena 10.44 -75.51 CO 6 0 .15
UIO Quito -0.13 -78.36 EC 5 0 .4
GYE Guayaquil -2.16 -79.88 EC 4 0 .4
LPB La_Paz -16.51 -68.19 BO 3 0 .2
VVI Santa_Cruz -17.64 -63.14 BO 3 0 .3
ASU Asuncion -25.24 -57.52 PY 1.5 0 1
MVD Montevideo -34.84 -56.03 UY 2.5 0 1
PTY Panama 9.07 -79.38 PA 18 .6 1
MEX Mexico_City 19.44 -99.07 MX 47 .1 .3
CUN Cancun 21.04 -86.87 MX 32 0 .7
GDL Guadalajara 20.52 -103.31 MX 18 0 .25
MTY Monterrey 25.78 -100.11 MX 14 0 .15
TIJ Tijuana 32.54 -116.97 MX 13 0 .05
SJD Los_Cabos 23.15 -109.72 MX 6 0 .6
MID Merida 20.94 -89.66 MX 4 0 .1
HMO Hermosillo 29.1 -111.05 MX 2.5 0 .05
"""

AU = """
SYD Sydney -33.95 151.18 AU 41 .1 .4
MEL Melbourne -37.67 144.84 AU 36 .05 .3
BNE Brisbane -27.38 153.12 AU 24 .05 .25
PER Perth -31.94 115.97 AU 15 0 .3
ADL Adelaide -34.95 138.53 AU 8.5 0 .1
OOL Gold_Coast -28.16 153.5 AU 6.5 0 .1
CBR Canberra -35.31 149.19 AU 3 0 0
CNS Cairns -16.88 145.75 AU 5 0 .1
DRW Darwin -12.41 130.88 AU 2 0 .1
HBA Hobart -42.84 147.51 AU 3 0 0
TSV Townsville -19.25 146.77 AU 1.8 0 0
ASP Alice_Springs -23.81 133.9 AU 0.6 0 0
AKL Auckland -37.01 174.79 NZ 20 .1 .5
"""

BLOC = {
    "GB": "EUR", "FR": "EUR", "DE": "EUR", "NL": "EUR", "BE": "EUR", "CH": "EUR", "AT": "EUR", "ES": "EUR",
    "PT": "EUR", "IT": "EUR", "GR": "EUR", "PL": "EUR", "CZ": "EUR", "HU": "EUR", "RO": "EUR", "BG": "EUR",
    "RS": "EUR", "DK": "EUR", "SE": "EUR", "NO": "EUR", "FI": "EUR", "IE": "EUR", "TR": "TR", "RU": "RU",
    "AZ": "CAS", "GE": "CAS", "AM": "CAS", "KZ": "CAS", "UZ": "CAS", "KG": "CAS", "AE": "GCC", "QA": "GCC",
    "BH": "GCC", "KW": "GCC", "SA": "GCC", "OM": "GCC", "JO": "LEV", "LB": "LEV", "IL": "LEV", "EG": "LEV",
    "IQ": "LEV", "IR": "LEV", "IN": "SAS", "PK": "SAS", "NP": "SAS", "BD": "SAS", "LK": "SAS", "TH": "SEA",
    "SG": "SEA", "MY": "SEA", "ID": "SEA", "PH": "SEA", "VN": "SEA", "MM": "SEA", "KH": "SEA", "LA": "SEA",
    "CN": "CHN", "HK": "HKT", "MO": "HKT", "TW": "HKT", "KR": "NEA", "JP": "NEA", "MN": "NEA",
    "ZA": "AFR", "KE": "AFR", "ET": "AFR", "NG": "AFR", "GH": "AFR", "CI": "AFR", "SN": "AFR", "TZ": "AFR",
    "UG": "AFR", "RW": "AFR", "AO": "AFR", "CD": "AFR", "ZM": "AFR", "ZW": "AFR", "MA": "NAF", "DZ": "NAF",
    "TN": "NAF", "CA": "CA", "BR": "SAM", "AR": "SAM", "CL": "SAM", "PE": "SAM", "CO": "SAM", "EC": "SAM",
    "BO": "SAM", "PY": "SAM", "UY": "SAM", "PA": "CAM", "MX": "CAM", "AU": "AU", "NZ": "NZ", "US": "US",
}

# Domestic calibration: national domestic segment pax/yr, enplanements per O&D trip, d0 (km), beta prior.
# Sources: G1-G10 (evidence.csv, unverified search summaries); others [K]/[U].
DOM = {
    "CN": (640e6, 1.12, 700, None, "G1"), "IN": (161e6, 1.12, 400, None, "G2"), "RU": (84e6, 1.25, 400, 0.5, "G3"),
    "BR": (93e6, 1.30, 400, None, "G4"), "AU": (63e6, 1.10, 400, None, "G10"), "CA": (None, 1.15, 400, 0.5, "G6"),
    "ZA": (18e6, 1.05, 400, 0.5, "G7+[K]"), "TR": (95e6, 1.15, 400, 0.5, "[K]"), "SA": (50e6, 1.1, 400, 0.5, "[K]"),
    "JP": (100e6, 1.1, 700, 0.5, "[K]"), "KR": (33e6, 1.0, 400, 0.5, "[K]"), "VN": (40e6, 1.05, 400, 0.5, "[K]"),
    "ID": (65e6, 1.15, 400, 0.5, "[K]"), "TH": (35e6, 1.1, 400, 0.5, "[K]"), "MY": (25e6, 1.1, 400, 0.5, "[K]"),
    "IR": (25e6, 1.05, 400, 0.5, "[K]"), "PK": (12e6, 1.05, 400, 0.5, "[K]"), "NG": (15e6, 1.05, 400, 0.5, "[K]"),
    "KE": (5e6, 1.1, 400, 0.5, "[K]"), "KZ": (10e6, 1.05, 400, 0.5, "[K]"), "EG": (4e6, 1.1, 400, 0.5, "[K]"),
    "GB": (20e6, 1.05, 500, 0.5, "[K]"), "FR": (25e6, 1.1, 600, 0.5, "[K]"), "DE": (20e6, 1.2, 600, 0.5, "[K]"),
    "ES": (45e6, 1.1, 500, 0.5, "[K]"), "IT": (35e6, 1.1, 500, 0.5, "[K]"), "CH": (1e6, 1.0, 400, 0.5, "[K]"),
    "AR": (13e6, 1.15, 400, 0.5, "[K]"), "CL": (15e6, 1.1, 400, 0.5, "[K]"), "PE": (13e6, 1.1, 400, 0.5, "[K]"),
    "CO": (28e6, 1.15, 400, 0.5, "[K]"), "EC": (3e6, 1.05, 400, 0.5, "[K]"), "BO": (3e6, 1.05, 400, 0.5, "[K]"),
    "MX": (60e6, 1.15, 400, 0.5, "[K]"), "PH": (25e6, 1.1, 400, 0.5, "[K]"), "NZ": (10e6, 1.1, 400, 0.5, "[K]"),
    "PE_": None,
}
COVER = 0.85      # share of national domestic O&D among the listed metros [U]

# Route anchors, O&D pax/yr = seats x LF 0.8 x local share (or pax / 1.0); metro codes.  src = evidence id.
ANCH_DOM = [
    ("CN", "PEK", "SHA", 7.71e6 * 0.8 * 0.95, "OAG24 via routes.csv"), ("CN", "CAN", "SHA", 7.01e6 * 0.8 * 0.95, "OAG24"),
    ("CN", "SHA", "SZX", 6.86e6 * 0.8 * 0.95, "OAG24"), ("CN", "PEK", "SZX", 5.60e6 * 0.8 * 0.95, "OAG24"),
    ("CN", "CAN", "HGH", 4.64e6 * 0.8 * 0.95, "OAG24"), ("CN", "CTU", "PEK", 4.60e6 * 0.8 * 0.95, "OAG24"),
    ("CN", "CAN", "PEK", 4.27e6 * 0.8 * 0.95, "OAG24"),
    ("IN", "BOM", "DEL", 7.96e6 * 0.8 * 0.9, "OAG24"), ("IN", "BLR", "DEL", 4.78e6 * 0.95, "FY24 pax"),
    ("IN", "BLR", "BOM", 4.43e6 * 0.95, "FY24 pax"), ("IN", "DEL", "HYD", 2.92e6 * 0.95, "FY24 pax"),
    ("IN", "CCU", "DEL", 2.82e6 * 0.95, "FY24 pax"),
    ("AU", "MEL", "SYD", 9.22e6 * 0.8 * 0.95, "OAG24"), ("AU", "BNE", "SYD", 4.6e6 * 0.95, "[K] BITRE pre-COVID"),
    ("AU", "BNE", "MEL", 3.6e6 * 0.95, "[K]"), ("AU", "MEL", "PER", 2.1e6 * 0.95, "[K]"),
    ("AU", "PER", "SYD", 1.9e6 * 0.95, "[K]"), ("AU", "ADL", "MEL", 2.4e6 * 0.95, "[K]"),
    ("CA", "YVR", "YYZ", 3.50e6 * 0.8 * 0.85, "G6"), ("CA", "YVR", "YYC", 2.96e6 * 0.8 * 0.8, "G6"),
    ("BR", "BSB", "SAO", 2.9e6 * 0.85, "G5"), ("BR", "CNF", "SAO", 2.3e6 * 0.9, "G5"), ("BR", "RIO", "SAO", 7.0e6, "G5 [U]"),
    ("ZA", "CPT", "JNB", 3.3e6, "G7"), ("SA", "JED", "RUH", 8.7e6 * 0.8 * 0.9, "OAG24"),
    ("VN", "HAN", "SGN", 10.63e6 * 0.8 * 0.95, "OAG24"), ("JP", "CTS", "TYO", 11.93e6 * 0.8, "OAG24"),
    ("JP", "FUK", "TYO", 11.34e6 * 0.8, "OAG24"),
]
ANCH_INT = [
    ("LON", "DXB", 2.70e6 * 0.8 * 0.6, "G9"), ("DXB", "RUH", 4.3e6 * 0.8 * 0.7, "G8"), ("CAI", "JED", 5.5e6 * 0.8 * 0.9, "G8"),
    ("BKK", "HKG", 4.2e6 * 0.8 * 0.85, "G8"), ("BKK", "SIN", 4.03e6 * 0.8 * 0.8, "G8"), ("CGK", "SIN", 4.07e6 * 0.8 * 0.8, "G8"),
    ("HKG", "TPE", 6.8e6 * 0.8 * 0.85, "G8"), ("SEL", "TYO", 5.0e6, "G8 + [U]"), ("OSA", "SEL", 4.98e6 * 0.8 * 0.85, "G8"),
    ("DXB", "BOM", 1.2e6, "[U]"), ("IST", "LON", 2.0e6, "[U]"), ("DEL", "LON", 0.9e6, "[U]"),
]
D0_INT = 400.0
# share of each bloc's international O&D that stays inside its modelled system [U]
INSYS = {"EUR": 0.75, "NEA": 0.6, "CHN": 0.8, "HKT": 0.75, "SEA": 0.85, "GCC": 0.9, "SAS": 0.9, "AFR": 0.9,
         "TR": 0.9, "RU": 0.95, "CAS": 0.95, "LEV": 0.95, "NAF": 0.95, "SAM": 0.45, "CAM": 0.3, "AU": 0.3, "NZ": 0.5}
# special-affinity bloc pairs [U]: Gulf-South Asia labour migration, Egypt/Levant-Saudi (pilgrimage, labour),
# Hong Kong-Taiwan.  Values set so the corresponding anchors land within ~x1.5 (calibration.csv).
AFF = {("GCC", "SAS"): 3.0, ("GCC", "LEV"): 3.0, ("HKT", "HKT"): 3.0}
# bloc-pair affinity groups calibrated on an anchor each; others = 1 (cross-checked by ANCH_INT residuals)
AFF_GROUPS = {("GCC", "SAS"): ["DXB-BOM"], ("GCC", "GCC"): ["DXB-RUH"], ("GCC", "LEV"): ["CAI-JED"],
              ("HKT", "SEA"): ["BKK-HKG"], ("SEA", "SEA"): ["BKK-SIN", "CGK-SIN"], ("HKT", "HKT"): ["HKG-TPE"],
              ("NEA", "NEA"): ["SEL-TYO", "OSA-SEL"], ("EUR", "TR"): ["IST-LON"], ("EUR", "SAS"): ["DEL-LON"],
              ("EUR", "GCC"): ["LON-DXB"]}


def parse(txt, system):
    rows = []
    for ln in txt.strip().splitlines():
        c, nm, la, lo, ct, pax, tr, ish = ln.split()
        rows.append(dict(code=c, name=nm.replace("_", " "), lat=float(la), lon=float(lo), country=ct,
                         pax_M=float(pax), transfer=float(tr), intl=float(ish), system=system))
    df = pd.DataFrame(rows)
    df["mass"] = df.pax_M * (1 - df.transfer) * 1e6
    df["m_dom"] = df.mass * (1 - df.intl)
    df["m_int"] = df.mass * df.intl
    df["bloc"] = df.country.map(BLOC)
    return df


def deter(d, beta, d0):
    d = np.maximum(d, 1.0)
    return d ** (-beta) * (1 - np.exp(-(d / d0) ** 2))


def gravity(df, log_rows):
    """Symmetric pair-trip matrix (both directions, trips/yr) for one system."""
    n = len(df)
    la, lo = df.lat.values, df.lon.values
    D = gc(la[:, None], lo[:, None], la[None, :], lo[None, :])
    idx = {c: i for i, c in enumerate(df.code)}
    T = np.zeros((n, n))
    ct = np.asarray(df.country.values, dtype=object)
    # --- domestic, per country
    for c in sorted(set(ct)):
        ii = np.where(ct == c)[0]
        if len(ii) < 2 or c not in DOM or DOM[c] is None:
            continue
        tot, ept, d0, bprior, src = DOM[c]
        md = df.m_dom.values[ii]
        Dc = D[np.ix_(ii, ii)]
        anch = [(idx[a], idx[b], v, s) for (cc, a, b, v, s) in ANCH_DOM if cc == c and a in idx and b in idx]
        best = None
        betas = [bprior] if bprior is not None else list(np.arange(-0.4, 1.61, 0.05))
        for beta in betas:
            base = np.outer(md, md) * deter(Dc, beta, d0)
            np.fill_diagonal(base, 0)
            res, w = [], []
            if tot:
                res.append(math.log(tot / ept * COVER / (base.sum() / 2)))
                w.append(2.0)
            for (a, b, v, s) in anch:
                ia, ib = np.where(ii == a)[0][0], np.where(ii == b)[0][0]
                res.append(math.log(v / base[ia, ib]))
                w.append(1.0)
            res, w = np.array(res), np.array(w)
            lk = (res * w).sum() / w.sum()
            err = float((w * (res - lk) ** 2).sum())
            if best is None or err < best[0]:
                best = (err, beta, math.exp(lk), base)
        err, beta, k, base = best
        Tc = k * base
        T[np.ix_(ii, ii)] += Tc
        for (a, b, v, s) in anch:
            ia, ib = np.where(ii == a)[0][0], np.where(ii == b)[0][0]
            log_rows.append(dict(kind="domestic", country=c, pair=f"{df.code[a]}-{df.code[b]}", km=round(D[a, b]),
                                 anchor_pax=v, model_pax=Tc[ia, ib], ratio=Tc[ia, ib] / v, source=s, beta=beta))
        log_rows.append(dict(kind="domestic total", country=c, pair="(all listed metros)", km=np.nan,
                             anchor_pax=(tot / ept * COVER) if tot else np.nan, model_pax=Tc.sum() / 2,
                             ratio=(Tc.sum() / 2 / (tot / ept * COVER)) if tot else np.nan, source=src, beta=beta))
    # --- international: production-constrained (symmetric IPF) gravity; each metro sends its in-system share
    #     of international O&D (m_int x INSYS[bloc]); beta fitted to anchors outside the affinity groups.
    mi = df.m_int.values * np.array([INSYS.get(b, 0.9) for b in df.bloc.values])
    bl = np.asarray(df.bloc.values, dtype=object)
    same_c = ct[:, None] == ct[None, :]
    aff = np.ones((n, n))
    for (g1, g2), f in AFF.items():
        m = ((bl[:, None] == g1) & (bl[None, :] == g2)) | ((bl[:, None] == g2) & (bl[None, :] == g1))
        aff[m] = f
    anch = [(idx[a], idx[b], v, s) for (a, b, v, s) in ANCH_INT if a in idx and b in idx]

    def ipf(beta):
        W = deter(D, beta, D0_INT) * aff
        W[same_c] = 0
        A = np.ones(n)
        for _ in range(200):
            row = (A[:, None] * A[None, :] * W * np.outer(mi, mi)).sum(1)
            A = A * np.sqrt(np.where(row > 0, mi / np.maximum(row, 1e-9), 1.0))
        return A[:, None] * A[None, :] * W * np.outer(mi, mi)

    best = None
    fit = [x for x in anch if not any({bl[x[0]], bl[x[1]]} == set(g) for g in AFF)]
    for beta in (np.arange(0.0, 2.01, 0.1) if fit else [1.0]):
        Ti = ipf(beta)
        err = sum(math.log(Ti[a, b] / v) ** 2 for a, b, v, s in fit) if fit else 0
        if best is None or err < best[0]:
            best = (err, beta, Ti)
    _, beta_i, Ti = best
    for (a, b, v, s) in anch:
        log_rows.append(dict(kind="international", country=f"{bl[a]}-{bl[b]}", pair=f"{df.code[a]}-{df.code[b]}",
                             km=round(D[a, b]), anchor_pax=v, model_pax=Ti[a, b], ratio=Ti[a, b] / v, source=s,
                             beta=beta_i))
    log_rows.append(dict(kind="international total", country="system", pair="sum", km=np.nan, anchor_pax=mi.sum() / 2,
                         model_pax=Ti.sum() / 2, ratio=Ti.sum() / mi.sum(), source="IPF", beta=beta_i))
    k = beta_i
    T = T + Ti
    np.fill_diagonal(T, 0)
    return D, T, k


def build_systems():
    rows = []
    sysd = {}
    ow = parse(OW, "OW")
    D, T, kint = gravity(ow, rows)
    globals()["_KINT"] = kint
    sysd["OW"] = dict(df=ow, D=D, T=T)
    for name, txt in (("CA", CA), ("LA", LA), ("AU", AU)):
        d = parse(txt, name)
        D, T, _ = gravity(d, rows)
        sysd[name] = dict(df=d, D=D, T=T)
    # US: DOT Table 1 + gravity fill (coverage / demand, read-only)
    inp = pickle.loads((NM / "demand" / ".cache" / "inputs.pkl").read_bytes())
    us = inp["df"].copy()
    us["country"] = "US"
    us["bloc"] = "US"
    us["system"] = "US"
    us["name"] = us.city
    la, lo = us.lat.values, us.lon.values
    D = gc(la[:, None], lo[:, None], la[None, :], lo[None, :])
    sysd["US"] = dict(df=us, D=D, T=np.asarray(inp["T"], float))
    cal = pd.DataFrame(rows)
    cal.to_csv(HERE / "calibration.csv", index=False, float_format="%.4g")
    return sysd


# =====================================================================================
# 2. Baselines
# =====================================================================================
SEA_WP = [(18, 64), (15, 68), (21, 65), (15, 88), (18, 89), (12, 86), (15, 113), (18, 116), (12, 110), (20, 108),
          (27, 51), (20, 39), (24, 36.5), (36, 18), (35, 25), (38, 5), (34, 30), (43, 34), (41, 51), (38.5, 25),
          (-10, -33), (-18, -37), (-24, -42), (-5, -36), (2, 3), (-5, 45), (-20, 42), (-20, -73), (-30, -73),
          (-10, -80), (13, -75), (25, -90), (27, -111), (-36, 130), (-32, 154), (-20, 155), (35, 123), (28, 124),
          (40, 134), (24, 119.5), (46, -8), (55, 3), (57, 19), (55, 148), (32, -75), (36, -125), (30, -79),
          (5, 95), (-37, 139), (8, 52), (12, 62), (25, -93), (22, -95), (16, -98)]


class Design:
    def __init__(self, name, M, a, b, cost_base, cost_int, fuel_ac, uplift_bl, uplift_sub):
        self.__dict__.update(locals())
        self.V = M * bt.A_MS * 3.6     # km/h


def water_time_min(Lw, D):
    """Minutes over an open-water run Lw (km) entered at M1.15; design Mach, kinematic ramps; 4 min FL600 extra."""
    v0 = bt.V_BL_START
    VM = D.M * bt.A_MS
    k = 1 / (2 * D.a) + 1 / (2 * D.b)
    A = np.maximum((Lw - bt.FL600_EXTRA[1]) * 1e3, 0)
    need = (VM ** 2 - v0 ** 2) * k
    Vp = np.where(A >= need, VM, np.sqrt(v0 ** 2 + A / k))
    cruise = np.maximum(A - (Vp ** 2 - v0 ** 2) * k, 0)
    t = ((Vp - v0) / D.a + (Vp - v0) / D.b + cruise / np.maximum(Vp, 1)) / 60 + bt.FL600_EXTRA[0]
    return np.where(Vp > v0 * 1.02, np.minimum(t, Lw / bt.V115 * 60), Lw / bt.V115 * 60)


def path_profile(pts, D, step=20.0):
    """pts: list of (lat, lon) waypoints. Returns (L_flown, land_cruise_km, water_minutes, restricted share)."""
    lats, lons, segs = [], [], []
    for p, q in zip(pts[:-1], pts[1:]):
        d = float(gc(p[0], p[1], q[0], q[1]))
        n = max(int(d / step), 1)
        f = (np.arange(n) + 0.5) / n
        la, lo = slerp(p, q, f)
        lats.append(la); lons.append(lo); segs.append(np.full(n, d / n))
    la, lo, ds = np.concatenate(lats), np.concatenate(lons), np.concatenate(segs) * bt.PATH
    s = np.cumsum(ds) - ds / 2
    L = s[-1] + ds[-1] / 2
    res = restricted(la, lo)
    a0 = bt.CLIMB[1] + bt.ACC_BL[1]
    a1 = L - bt.APPR[1] - bt.DESC_SUB[1] - bt.DEC_BL[1]
    win = (s >= a0) & (s <= a1)
    land_km, water_min = 0.0, 0.0
    if win.any():
        r, d = res[win], ds[win]
        # runs
        brk = np.flatnonzero(np.diff(r.astype(int))) + 1
        for rr, dd in zip(np.split(r, brk), np.split(d, brk)):
            Lr = dd.sum()
            if rr[0]:
                land_km += Lr
            else:
                tw = float(water_time_min(np.array([Lr]), D)[0])
                if tw < Lr / bt.V115 * 60 - 0.01:
                    water_min += tw
                else:
                    land_km += Lr
    return L, land_km, water_min, float((res * ds).sum() / L)


def t_bl_from(L, land_km, water_min, V_bl_kmh=None):
    V = V_bl_kmh or bt.V115
    t = np.where(L >= 695, bt.TAXI + 57 + land_km / V * 60 + water_min, np.inf)
    sub = bt.t_direct(L / bt.PATH, "sub")
    o94 = np.where(L >= 490, bt.TAXI + 46 + (L - 490) / bt.V094 * 60, sub)
    return np.minimum(np.minimum(t, o94), sub)


def baselines(S, D):
    df, Dm, T = S["df"], S["D"], S["T"]
    iu, ju = np.triu_indices(len(df), 1)
    d, pax = Dm[iu, ju], T[iu, ju]
    keep = (d >= 1000) & (pax >= 2000)
    iu, ju, d, pax = iu[keep], ju[keep], d[keep], pax[keep]
    la, lo = df.lat.values, df.lon.values
    wp = np.array(SEA_WP, float)
    out = dict(iu=iu, ju=ju, d=d, pax=pax, L=np.zeros(len(iu)), land=np.zeros(len(iu)), water=np.zeros(len(iu)),
               rfrac=np.zeros(len(iu)), via=np.full(len(iu), -1))
    t0 = time.time()
    for n, (i, j) in enumerate(zip(iu, ju)):
        p, q = (la[i], lo[i]), (la[j], lo[j])
        best = path_profile([p, q], D)
        tb = float(t_bl_from(np.array(best[0]), best[1], best[2]))
        via = -1
        if d[n] >= 1300:
            dd = gc(p[0], p[1], wp[:, 0], wp[:, 1]) + gc(wp[:, 0], wp[:, 1], q[0], q[1])
            for w in np.where(dd <= 1.2 * d[n])[0]:
                pr = path_profile([p, tuple(wp[w]), q], D)
                t = float(t_bl_from(np.array(pr[0]), pr[1], pr[2]))
                if t < tb - 0.5:
                    tb, best, via = t, (pr[0], pr[1], pr[2], best[3]), w
        out["L"][n], out["land"][n], out["water"][n], out["rfrac"][n] = best[0], best[1], best[2], best[3]
        out["via"][n] = via
        if n % 2000 == 0:
            log(f"   baseline {n}/{len(iu)}  {time.time() - t0:.0f}s")
    out["t_bl"] = t_bl_from(out["L"], out["land"], out["water"])
    out["t_sub"] = bt.t_direct(d, "sub")
    ub = np.minimum(*(bt.t_trunk(0.0, d, 0.0, D.M, f, 5.0, D.a, D.b)[0] for f in ("sub", "bl")))
    out["t_own"] = ub
    return out


# =====================================================================================
# 3. Trunks and routing
# =====================================================================================
class Trunk:
    def __init__(self, edges, P, step=25.0):
        """edges: list of (name_a, name_b); P: dict name -> (lat, lon)."""
        G = nx.Graph()
        pts, key = [], {}

        def pid(nm):
            if nm not in key:
                key[nm] = len(pts)
                pts.append(P[nm])
            return key[nm]
        self.edge_pts = []
        for a, b in edges:
            ia, ib = pid(a), pid(b)
            pa, pb = P[a], P[b]
            dkm = float(gc(pa[0], pa[1], pb[0], pb[1]))
            n = max(int(round(dkm / step)), 1)
            prev = ia
            chain = [ia]
            if n > 1:
                la, lo = slerp(pa, pb, np.arange(1, n) / n)
                for x, y in zip(la, lo):
                    k = len(pts)
                    pts.append((float(x), float(y)))
                    G.add_edge(prev, k, length=dkm / n)
                    prev = k
                    chain.append(k)
            G.add_edge(prev, ib, length=dkm / n)
            chain.append(ib)
            self.edge_pts.append((a, b, chain))
        self.G = G
        self.ll = np.array(pts)
        n = len(pts)
        Dt = np.full((n, n), np.inf)
        for s, dist in nx.all_pairs_dijkstra_path_length(G, weight="length"):
            for t, v in dist.items():
                Dt[s, t] = v
        self.Dt = Dt
        # segment lengths and restricted flags (midpoints)
        u = np.array([e[0] for e in G.edges]); v = np.array([e[1] for e in G.edges])
        self.e_u, self.e_v = u, v
        self.e_len = np.array([G[a][b]["length"] for a, b in zip(u, v)])
        mid = (self.ll[u] + self.ll[v]) / 2
        self.e_res = restricted(mid[:, 0], mid[:, 1])
        self.km = float(self.e_len.sum())
        self.km_res = float((self.e_len * self.e_res).sum())
        self.edges = list(edges)

    def ramp_dist(self, la, lo, R=16):
        Dm = gc(la[:, None], lo[:, None], self.ll[None, :, 0], self.ll[None, :, 1])
        R = min(R, Dm.shape[1])
        idx = np.argsort(Dm, axis=1)[:, :R]
        return idx, np.take_along_axis(Dm, idx, 1)


def route(tr, df, B, D, sel=None, R=16):
    """Best same-aircraft time via the trunk for pairs in B (optionally a subset)."""
    sel = np.arange(len(B["iu"])) if sel is None else sel
    la, lo = df.lat.values, df.lon.values
    idx, dist = tr.ramp_dist(la, lo, R)
    iu, ju = B["iu"][sel], B["ju"][sel]
    t_best = np.full(len(sel), np.inf)
    on = np.zeros(len(sel), int); off = np.zeros(len(sel), int)
    CH = 4000
    for s in range(0, len(sel), CH):
        a, b = iu[s:s + CH], ju[s:s + CH]
        ai, ad = idx[a][:, :, None], dist[a][:, :, None]
        bi, bd = idx[b][:, None, :], dist[b][:, None, :]
        L = tr.Dt[ai, bi]
        best = np.full(L.shape, np.inf)
        for f in ("sub", "bl"):
            t, _ = bt.t_trunk(ad, L, bd, D.M, f, 5.0, D.a, D.b)
            best = np.minimum(best, t)
        flat = best.reshape(len(a), -1)
        k = np.argmin(flat, 1)
        t_best[s:s + CH] = flat[np.arange(len(a)), k]
        Rr = idx.shape[1]
        on[s:s + CH] = idx[a, k // Rr]
        off[s:s + CH] = idx[b, k % Rr]
    return t_best, on, off


def evaluate(tr, S, B, D, R=16, loads=False, V_bl=None, rng_km=None, mask=None):
    ok = (B["t_sub"] - B["t_own"] >= CUT_MIN) & (B["d"] <= (rng_km or RANGE_KM))
    if mask is not None:
        ok &= mask
    sel = np.where(ok)[0]
    t, on, off = route(tr, S["df"], B, D, sel, R)
    t_bl = B["t_bl"][sel] if V_bl is None else t_bl_from(B["L"][sel], B["land"][sel], B["water"][sel], V_bl)
    pax = B["pax"][sel]
    s_bl = t_bl - t
    s_bl = np.where(s_bl >= CUT_MIN, s_bl, 0.0)
    s_sub = B["t_sub"][sel] - t
    s_sub = np.where(s_sub >= CUT_MIN, s_sub, 0.0)
    r = dict(sel=sel, t=t, pax=pax, s_bl=s_bl, s_sub=s_sub, on=on, off=off,
             B_bl=float((pax * s_bl).sum() / 60), B_sub=float((pax * s_sub).sum() / 60),
             users_bl=float(pax[s_bl > 0].sum()), users_sub=float(pax[s_sub > 0].sum()),
             km=tr.km, km_res=tr.km_res)
    if loads:
        r["load"] = edge_loads(tr, on, off, pax, s_bl > 0)
        r["load_sub"] = edge_loads(tr, on, off, pax, s_sub > 0)
    return r


def edge_loads(tr, on, off, pax, use):
    load = {}
    cache = {}
    agg = {}
    for a, b, p in zip(on[use], off[use], pax[use]):
        if a == b:
            continue
        kk = (min(a, b), max(a, b))
        agg[kk] = agg.get(kk, 0.0) + p
    for kk, p in agg.items():
        path = nx.shortest_path(tr.G, kk[0], kk[1], weight="length")
        for u, v in zip(path[:-1], path[1:]):
            e = (min(u, v), max(u, v))
            load[e] = load.get(e, 0.0) + p
    return load


# =====================================================================================
# 4. Corridors (hand-built spines).  Waypoints W_* are inland Steiner points [author choice].
# =====================================================================================
WAYPOINTS = {
    "W_ANATOLIA": (39.0, 36.0), "W_IRAQ": (33.0, 44.5), "W_IRAN_S": (28.5, 56.5), "W_BALUCH": (27.0, 63.0),
    "W_SINDH": (26.5, 69.0), "W_BURMA": (22.0, 95.5), "W_LAOS": (20.5, 102.0), "W_GUANGXI": (23.5, 107.5),
    "W_ODISHA": (21.0, 84.0), "W_MP": (23.3, 77.4), "W_KAZ_W": (51.0, 58.0), "W_KAZ_E": (47.0, 76.0),
    "W_GANSU": (39.5, 98.5), "W_MONG": (45.5, 108.0), "W_TURKMEN": (38.0, 58.4), "W_AFGH": (34.5, 66.0),
    "W_NULLARBOR": (-31.0, 129.5), "W_SA_N": (-31.0, 136.5), "W_SUDAN": (15.6, 32.5), "W_TANZ": (-6.2, 35.7),
    "W_ZAMBIA": (-13.0, 28.6), "W_CONGO": (-4.4, 15.3), "W_CAMEROON": (3.9, 11.5), "W_BAHIA": (-12.5, -41.0),
    "W_PIAUI": (-7.5, -42.0), "W_TOCANT": (-10.0, -48.3), "W_RONDON": (-12.0, -61.0), "W_ONT_N": (49.5, -86.0),
    "W_SASK": (50.5, -105.0), "W_BC": (50.5, -119.0), "W_VOLGA": (55.5, 47.0), "W_CHITA": (52.0, 113.5),
    "W_AMUR": (51.0, 126.0), "W_CHACO": (-26.0, -60.0), "W_ANDES": (-32.7, -69.5), "W_SAHEL": (12.0, 14.0),
    "W_BALKAN": (44.0, 21.0), "W_JORDAN": (30.5, 38.0), "W_NAJD": (27.0, 44.0), "W_HUBEI": (31.0, 112.0),
    "W_HENAN": (34.0, 114.0), "W_SICHUAN_N": (32.5, 105.0), "W_QINGHAI": (37.0, 100.0), "W_HAMI": (42.8, 93.5),
    "W_GUIZHOU": (26.5, 106.7), "W_DECCAN": (17.5, 76.5), "W_UP": (26.0, 82.0), "W_CHHATT": (21.5, 81.6),
    "W_KANSAS": (38.8, -98.0), "W_SONORA": (29.0, -110.0), "W_DURANGO": (24.5, -104.5),
}

# name, region, edges (chains "A-B-C" or explicit pairs), note
CORRIDORS = [
    ("US north spine (reference)", "US", ["LAS-DEN-OMA-ORD-CLE-LGA"], "coverage.md spine"),
    ("US southern spine", "US", ["LAX-PHX-DFW-ATL-DCA"], "coverage.md 'south only'"),
    ("Canada YVR-YYC-YWG-YYZ-YUL", "CA", ["YVR-YYC-YWG-YYZ-YUL"], "GC YWG-YYZ crosses Lake Superior"),
    ("Canada YVR-YYC-YYZ-YUL (direct)", "CA", ["YVR-YYC-YYZ-YUL"], ""),
    ("Australia PER-ADL-MEL-SYD", "AU", ["PER-ADL-MEL-SYD"], "PER-ADL crosses the Bight (free water)"),
    ("Australia PER-SYD + MEL, BNE arms", "AU", ["PER-W_SA_N-SYD", "W_SA_N-MEL", "SYD-BNE"], ""),
    ("Australia PER-SYD-BNE", "AU", ["PER-SYD-BNE"], ""),
    ("LON-IST-DXB-BOM-HKG (user)", "OW", ["LON-W_BALKAN-IST-W_ANATOLIA-W_IRAQ-DXB-BOM-W_ODISHA-W_BURMA-W_LAOS-HKG"],
     "Gulf + Arabian Sea + Bay of Bengal gaps"),
    ("FRA-IST-DXB-DEL-CCU-BKK", "OW", ["LON-FRA-W_BALKAN-IST-W_ANATOLIA-W_IRAQ-DXB-W_IRAN_S-W_BALUCH-W_SINDH-DEL-CCU-W_BURMA-BKK"],
     "overland through Iran/Pakistan"),
    ("LON-IST-DXB, Europe-Gulf only", "OW", ["LON-FRA-W_BALKAN-IST-W_ANATOLIA-W_IRAQ-DXB"], ""),
    ("Europe-Gulf-India (LON-IST-DXB-DEL)", "OW", ["LON-FRA-W_BALKAN-IST-W_ANATOLIA-W_IRAQ-DXB-W_IRAN_S-W_BALUCH-W_SINDH-DEL"], ""),
    ("China N-S PEK-WUH-CAN-SZX/HKG", "OW", ["PEK-W_HENAN-WUH-CAN-SZX", "CAN-HKG"], ""),
    ("China E-W SHA-WUH-CKG-CTU", "OW", ["SHA-WUH-CKG-CTU"], ""),
    ("China X: PEK-CAN + SHA-CTU", "OW", ["PEK-W_HENAN-WUH-CAN-SZX", "CAN-HKG", "SHA-WUH-CKG-CTU"], "crossing at Wuhan"),
    ("China X + CTU-KMG", "OW", ["PEK-W_HENAN-WUH-CAN-SZX", "CAN-HKG", "SHA-WUH-CKG-CTU-KMG"], ""),
    ("China NW: PEK-LHW-URC + SHA-XIY-LHW", "OW", ["PEK-W_HENAN-XIY-LHW-W_GANSU-W_HAMI-URC", "SHA-XIY"], "Gansu corridor"),
    ("China big H: PEK-CAN, SHA-CTU-KMG, CTU-URC", "OW", ["PEK-W_HENAN-WUH-CAN-SZX", "CAN-HKG", "SHA-WUH-CKG-CTU-KMG",
                                                      "CTU-W_SICHUAN_N-LHW-W_GANSU-W_HAMI-URC"], ""),
    ("India DEL-BOM-BLR", "OW", ["DEL-BOM-BLR"], ""),
    ("India DEL-HYD-BLR/MAA + CCU", "OW", ["CCU-W_UP-DEL", "DEL-W_MP-HYD-BLR", "HYD-MAA"], ""),
    ("India Y: DEL-MP-{BOM, HYD-BLR} + CCU", "OW", ["CCU-W_UP-DEL", "DEL-W_MP", "W_MP-BOM", "W_MP-HYD-BLR"], ""),
    ("Trans-Siberian LON-MOW-OVB-PEK/SEL", "OW", ["LON-AMS-BER-WAW-MOW-W_VOLGA-SVX-OVB-KJA-IKT-W_MONG-PEK-SEL"],
     "Russian airspace closed to EU/US carriers"),
    ("Europe-China via Kazakhstan", "OW", ["LON-FRA-WAW-MOW-W_KAZ_W-W_KAZ_E-URC-W_HAMI-W_GANSU-LHW-XIY-PEK"], ""),
    ("Gulf-China overland DXB-URC-PEK", "OW", ["DXB-W_IRAN_S-W_TURKMEN-TAS-ALA-URC-W_HAMI-W_GANSU-LHW-XIY-PEK"],
     "Iran, Central Asia"),
    ("Gulf-India-China DXB-DEL-KMG-CAN", "OW", ["DXB-W_IRAN_S-W_BALUCH-W_SINDH-DEL-W_UP-KTM-W_BURMA-KMG-W_GUANGXI-CAN-HKG"],
     "Himalayan foothills, Myanmar"),
    ("Russia MOW-SVX-OVB-IKT-KHV-VVO", "OW", ["MOW-W_VOLGA-SVX-OVB-KJA-IKT-W_CHITA-W_AMUR-KHV-VVO"], ""),
    ("Russia MOW-SVX-OVB", "OW", ["MOW-W_VOLGA-SVX-OVB"], ""),
    ("China-SE Asia KMG-BKK-KUL-SIN", "OW", ["CTU-KMG-W_LAOS-BKK-KUL-SIN"], "Malay peninsula coastal"),
    ("China-SE Asia PEK-CAN-HAN-BKK-SIN", "OW", ["PEK-W_HENAN-WUH-CAN-W_GUANGXI-HAN-W_LAOS-BKK-KUL-SIN"], ""),
    ("Africa CAI-KRT-ADD-NBO-JNB", "OW", ["CAI-W_SUDAN-ADD-NBO-W_TANZ-W_ZAMBIA-JNB"], "Sudan war"),
    ("Africa LOS-FIH-JNB", "OW", ["LOS-W_CAMEROON-W_CONGO-LAD-W_ZAMBIA-JNB"], ""),
    ("Brazil SAO-BSB-FOR/REC", "LA", ["SAO-BSB-W_PIAUI-FOR", "W_PIAUI-REC"], ""),
    ("Brazil SAO-BSB-MAO + BEL", "LA", ["SAO-BSB-W_RONDON-MAO", "BSB-W_TOCANT-BEL"], "Amazon"),
    ("Brazil SAO-BSB + NE + Amazon", "LA", ["RIO-SAO-BSB-W_PIAUI-FOR", "W_PIAUI-REC", "BSB-W_TOCANT-BEL", "BSB-W_RONDON-MAO"], ""),
    ("S. Cone SCL-BUE-SAO (Andes)", "LA", ["SCL-W_ANDES-COR-BUE-SAO"], "Andes crossing"),
    ("Andes BOG-LIM-SCL (inland)", "LA", ["BOG-UIO-LIM-LPB-SCL"], "Pacific offshore alternative"),
    ("Mexico TIJ-GDL-MEX-CUN", "LA", ["TIJ-W_SONORA-GDL-MEX-CUN"], ""),
]


def chains_to_edges(chains):
    out = []
    for c in chains:
        nodes = c.split("-")
        # tokens may be W_ names containing '_' only
        out += list(zip(nodes[:-1], nodes[1:]))
    return out


# =====================================================================================
# 5. Greedy builder (trunk.md style)
# =====================================================================================
SCOPES = {
    "US": ("US", None), "Canada": ("CA", None), "Australia": ("AU", ["AU"]), "Brazil": ("LA", ["BR"]),
    "China": ("OW", ["CN", "HK"]), "India": ("OW", ["IN"]), "Russia": ("OW", ["RU"]),
    "Europe-Middle East-South Asia": ("OW", ["EUR", "TR", "GCC", "LEV", "SAS"]),
    "Africa": ("OW", ["AFR", "LEV"]), "China-SE Asia": ("OW", ["CN", "HK", "SEA"]),
}


def greedy(S, B, D, P, members, theta=0.5, max_steps=7, n_seed=3, logf=None):
    df = S["df"]
    nodes = [c for c in df.code[members]]
    # cut tiny nodes for speed: keep top 45 by mass
    mass = df["mass"].values if "mass" in df else S["T"].sum(1)
    m = pd.Series(mass, index=df.code.values).loc[nodes].sort_values(ascending=False)
    nodes = list(m.index[:45])
    ll = np.array([P[c] for c in nodes])
    xy = np.c_[ll[:, 1] * np.cos(np.radians(ll[:, 0].mean())), ll[:, 0]]
    tri = Delaunay(xy)
    G = nx.Graph()
    for s in tri.simplices:
        for i in range(3):
            a, b = nodes[s[i]], nodes[s[(i + 1) % 3]]
            d = float(gc(P[a][0], P[a][1], P[b][0], P[b][1]))
            if d <= 1800:
                G.add_edge(a, b, length=d)
    # seeds: densest pairs by own-GC potential within scope
    code = np.asarray(df.code.values, dtype=object)
    inn = np.isin(code[B["iu"]], nodes) & np.isin(code[B["ju"]], nodes)
    pot = B["pax"] * np.maximum(B["t_bl"] - B["t_own"], 0) * inn * (B["d"] <= RANGE_KM)
    seeds = []
    for k in np.argsort(-pot)[:20]:
        a, b = code[B["iu"][k]], code[B["ju"][k]]
        if a in G and b in G and nx.has_path(G, a, b):
            pth = nx.shortest_path(G, a, b, weight="length")
            e = list(zip(pth[:-1], pth[1:]))
            if e not in seeds:
                seeds.append(e)
        if len(seeds) >= n_seed:
            break
    best_overall = None
    allm = set(df.code[members])
    pm = np.isin(code[B["iu"]], list(allm)) & np.isin(code[B["ju"]], list(allm))
    for seed in seeds:
        edges = list(seed)
        tr = Trunk(edges, P)
        r = evaluate(tr, S, B, D, R=10, mask=pm)
        val = r["B_bl"] / max(tr.km_res, 50)
        hist = [(list(edges), r["B_bl"], tr.km_res)]
        for step in range(max_steps):
            tn = {a for e in edges for a in e}
            dist, paths = nx.multi_source_dijkstra(G, tn, weight="length")
            cands = [(n_, paths[n_]) for n_ in dist if n_ not in tn and dist[n_] <= 2200]
            cands = sorted(cands, key=lambda x: -m[x[0]])[:18]
            best = None
            for n_, pth in cands:
                e2 = edges + list(zip(pth[:-1], pth[1:]))
                t2 = Trunk(e2, P)
                r2 = evaluate(t2, S, B, D, R=10, mask=pm)
                dk = max(t2.km_res - tr.km_res, 25.0)
                marg = (r2["B_bl"] - r["B_bl"]) / dk
                if best is None or marg > best[0]:
                    best = (marg, e2, t2, r2)
            if best is None:
                break
            avg = r["B_bl"] / max(tr.km_res, 50)
            if best[0] < theta * avg:
                break
            _, edges, tr, r = best
            hist.append((list(edges), r["B_bl"], tr.km_res))
            if logf:
                logf(f"    step {step}: {len(edges)} edges, {tr.km_res:.0f} km_res, B_bl {r['B_bl'] / 1e6:.2f} M, marg {best[0]:.0f} vs avg {avg:.0f}")
        val = r["B_bl"] / max(tr.km_res, 50)
        if best_overall is None or val > best_overall[0]:
            best_overall = (val, edges, hist)
    return best_overall


# =====================================================================================
# 6. Climate: Mach cutoff and jet exposure (parametric climatology [U])
# =====================================================================================
def t_sfc(lat, lon, month, elev):
    """Monthly mean surface temperature, K [U]: zonal fit + continental amplitude + hot-desert bonus - lapse."""
    alat = np.abs(lat)
    Tann = 28 - 0.0076 * alat ** 2
    cont = np.where((lon > -12) & (lon < 25) & (lat > 40) & (lat < 62), 0.6, 1.0)         # maritime Europe
    cont = np.where((lat < 0), 0.75, cont)                                               # southern hemisphere
    amp = 0.35 * alat * cont
    desert = ((lat > 15) & (lat < 33) & (lon > -15) & (lon < 62)) | ((lat < -20) & (lat > -32) & (lon > 115) & (lon < 145))
    Tann = Tann + np.where(desert, 4.0, 0.0)
    ph = np.where(lat >= 0, np.cos(2 * np.pi * (month - 7) / 12), np.cos(2 * np.pi * (month - 1) / 12))
    return 273.15 + Tann + amp * ph - 6.5 * (elev - 0.3)


def jet_u(lat, lon, month):
    """Zonal wind at ~11.5 km, m/s [U]: subtropical + polar-front jets, tropical easterly jet (JJA)."""
    nh = lat >= 0
    w = np.where(nh, (1 + np.cos(2 * np.pi * (month - 1) / 12)) / 2, (1 + np.cos(2 * np.pi * (month - 7) / 12)) / 2)
    alat = np.abs(lat)
    lon_ = ((lon + 180) % 360) - 180
    amp = np.full(np.shape(lat), 30.0)
    amp = np.where(nh & (lon_ >= 100) & (lon_ <= 150), 60.0, amp)
    amp = np.where(nh & (lon_ >= 30) & (lon_ < 100), 45.0, amp)
    amp = np.where(nh & (lon_ >= -130) & (lon_ <= -60), 38.0, amp)
    amp = np.where(~nh & (((lon_ >= 110) & (lon_ <= 160)) | ((lon_ >= -80) & (lon_ <= -40))), 45.0, amp)
    phi_s = 30 + 10 * (1 - w)
    U_s = amp * (0.4 + 0.6 * w) * np.exp(-((alat - phi_s) / 12.0) ** 2)
    U_p = 25.0 * np.exp(-((alat - (50 + 5 * (1 - w))) / 12.0) ** 2)
    tej = np.where(nh & (lon_ > 30) & (lon_ < 110), -25.0 * (1 - w) * np.exp(-((lat - 13) / 6.0) ** 2), 0.0)
    return U_s + U_p, tej


def cutoff_stats(pts_ll, brg, rng, n=600):
    """Per point: draw (month, direction, weather); Mach cutoff at FL470 (14.3 km).  Returns dict."""
    lat, lon = pts_ll[:, 0], pts_ll[:, 1]
    elev = elevation(lat, lon)
    month = rng.integers(1, 13, (len(lat), n))
    dirn = rng.choice([1.0, -1.0], (len(lat), n))
    LAT, LON, EL = lat[:, None], lon[:, None], elev[:, None]
    Ts = t_sfc(LAT, LON, month, EL) + rng.normal(0, 1, (len(lat), n)) * np.where(np.abs(LAT) < 20, 2.0, 5.0)
    uj, tej = jet_u(LAT, LON, month)
    uj = uj + rng.normal(0, 12.0, uj.shape)
    alat = np.abs(LAT)
    Th = 211 + 6 * np.minimum(alat, 30) / 30        # K at 14.3 km
    h = 14.3
    comp = np.sin(np.radians(brg))[:, None] * dirn   # eastward component of the track direction
    z = np.linspace(0, h, 30)[None, None, :]
    zz = EL[..., None] + (h - EL[..., None]) * z / h
    T = np.maximum(Ts[..., None] - 6.5 * (zz - EL[..., None]), Th[..., None])
    u = (uj[..., None] * np.exp(-((zz - 11.5) / 4.0) ** 2) + tej[..., None] * np.exp(-((zz - 15.0) / 3.0) ** 2)
         + 5.0 * np.exp(-zz / 2.0))
    a = 20.05 * np.sqrt(T)
    tail = comp[..., None] * u
    a_h = 20.05 * np.sqrt(Th)
    u_h = comp * (uj * math.exp(-((h - 11.5) / 4.0) ** 2) + tej * math.exp(-((h - 15) / 3.0) ** 2))
    Mc = ((a + tail).max(-1) - u_h) / a_h
    return dict(Mc=Mc, share115=float((Mc >= 1.15).mean()), mean_eff=float(np.clip(Mc - 0.02, 0.94, 1.3).mean()),
                Mc_median=float(np.median(Mc)), jet_winter_ms=float(np.abs(jet_u(lat, lon, np.where(lat >= 0, 1, 7))[0]).mean()),
                jet_u2=float((jet_u(LAT, LON, month)[0] ** 2).mean()), elev_max=float(elev.max()),
                elev_mean=float(elev.mean()))


# =====================================================================================
# 7. Economics
# =====================================================================================
def designs():
    md = pd.read_csv(NM / "highmach" / "mach_designs.csv").set_index("design")
    I = pd.read_csv(BS / "economics" / "integration" / "results.csv")
    m = I[I.design.str.startswith("M2.0")].set_index("scenario").net_musd_km_yr
    int_ratio = m["integrated, aggressive levers, central"] / m["cost-down aggressive + 2nd tether row"]
    bc = pd.read_csv(NM / "highmach" / "bcr_by_mach.csv")
    out = {}
    for key, nm, eng in (("M2.2", "A: M2.2 Al-limit", "adaptive"), ("M2.0", "M2.0 adaptive", "adaptive")):
        r = md.loc[nm]
        s = bc[(bc.design == nm) & (bc.network == "spine") & (bc.demand == "base") & (bc.capture == 1.0) & (bc.cost == "base")]
        up = {vs: float(s[s.vs == vs].bcr.iloc[0] / s[s.vs == vs].bcr_time.iloc[0]) for vs in ("bl", "sub")}
        fuel_ac = ha.size(float(r.M), eng, kT=1.40)
        out[key] = Design(key, float(r.M), float(r.a_acc), float(r.a_dec), float(r.cost_base_musd_km_yr),
                          float(r.cost_aggressive_musd_km_yr) * int_ratio, fuel_ac, up["bl"], up["sub"])
    ref = ha.size(1.7, "fixed", kT=1.22)
    return out, ref, int_ratio


def fuel_per_pax(ac, km):
    L = np.asarray(km) * 1.03
    fa = ac["f_accel"]
    f = ac["W0_t"] * 1e3 * (fa + (1 - fa) * (1 - np.exp(-L / ac["RF_km"]))) * 0.97
    return f / (ha.SEATS * LF)


def econ(r, B, D, ref, vot_factor=1.0):
    """B/C rows at capture 1 (scale linearly).  Gross includes highmach's induced/day-return uplift."""
    sel = r["sel"]
    d = B["d"][sel]
    pax = r["pax"]
    usd_per_kg = (FUEL_USD_KG + CO2_KG / 1e3 * SCC) * AVAIL
    use_bl, use_sub = r["s_bl"] > 0, r["s_sub"] > 0
    f_M = fuel_per_pax(D.fuel_ac, d)
    pen_bl = float((pax * use_bl * (f_M - fuel_per_pax(ref, d))).sum() * usd_per_kg)
    pen_sub = float((pax * use_sub * (f_M - FUEL_SUB_PAX_KM * d * 1.03)).sum() * usd_per_kg)
    ben_bl = r["B_bl"] * VOT * vot_factor * AVAIL * D.uplift_bl
    ben_sub = r["B_sub"] * VOT * vot_factor * AVAIL * D.uplift_sub
    out = {}
    for lvl, c in (("base", D.cost_base), ("int", D.cost_int)):
        cost = max(r["km_res"], 1.0) * c * 1e6
        out[f"bc_bl_{lvl}"] = ben_bl / cost
        out[f"bc_bl_{lvl}_net"] = (ben_bl - pen_bl) / cost
        out[f"bc_sub_{lvl}"] = ben_sub / cost
        out[f"bc_sub_{lvl}_net"] = (ben_sub - pen_sub) / cost
    out["fuel_co2_pen_bl_busd"] = pen_bl / 1e9
    return out


GDPPC = {  # GDP per capita PPP 2024, $k [K, approximate]
    "US": 86, "CA": 62, "AU": 68, "NZ": 52, "GB": 58, "FR": 60, "DE": 70, "NL": 75, "BE": 70, "CH": 90, "AT": 72,
    "ES": 52, "PT": 45, "IT": 58, "GR": 40, "PL": 50, "CZ": 52, "HU": 45, "RO": 45, "BG": 35, "RS": 27, "DK": 80,
    "SE": 70, "NO": 90, "FI": 63, "IE": 90, "TR": 42, "RU": 45, "AZ": 23, "GE": 25, "AM": 20, "KZ": 37, "UZ": 11,
    "KG": 7, "AE": 80, "QA": 120, "BH": 60, "KW": 55, "SA": 70, "OM": 40, "JO": 11, "LB": 15, "IL": 55, "EG": 18,
    "IQ": 12, "IR": 19, "IN": 10.5, "PK": 6.5, "NP": 5, "BD": 9.6, "LK": 15, "TH": 23, "SG": 140, "MY": 37,
    "ID": 16, "PH": 11, "VN": 15, "MM": 5.5, "KH": 6, "LA": 10, "CN": 26, "HK": 75, "MO": 80, "TW": 75, "KR": 57,
    "JP": 52, "MN": 18, "ZA": 16, "KE": 6.5, "ET": 3.5, "NG": 6, "GH": 7, "CI": 7, "SN": 4.5, "TZ": 4, "UG": 3,
    "RW": 3.2, "AO": 8, "CD": 1.5, "ZM": 4, "ZW": 3, "MA": 10, "DZ": 16, "TN": 13, "BR": 22, "AR": 30, "CL": 33,
    "PE": 17, "CO": 21, "EC": 15, "BO": 11, "PY": 17, "UY": 32, "PA": 40, "MX": 25,
}


# =====================================================================================
# main
# =====================================================================================
FLAGS = {  # airspace / politics flag (text only) and jurisdiction for regulation
    "US": ("US/FAA: §91.817 ban; 0.11 psf NPRM (Jul 2026)", "low"),
    "CA": ("Canada CARs: Mach >= 1 banned everywhere [R2]", "low"),
    "AU": ("Australia: no rule text found; Concorde barred overland [R5]", "low"),
    "CN": ("China CAAC: rule not found [R6]; military-controlled airspace (~20-30 % civil)", "medium"),
    "IN": ("India DGCA: rule not found [R6]; IAF-shared airspace", "low-medium"),
    "RU": ("Russia: Tu-144 flew overland 1977-78 [K1]; sanctions, airspace closed to EU/US/UK carriers", "high"),
    "BR": ("Brazil ANAC/DECEA: rule not found", "low"),
    "EU": ("UK CAA proposing overland ban [R3]; EASA silent; national bans likely [R4]", "low (EU) / medium"),
}


def main(stage_max=9):
    global GRID
    GRID = restricted_grid()
    log("grid ok", GRID.shape, GRID.mean())
    DS, REF, int_ratio = designs()
    log("designs", {k: (v.M, round(v.a, 3), round(v.cost_base, 2), round(v.cost_int, 2), round(v.uplift_bl, 3)) for k, v in DS.items()})
    SYS = cached("systems", build_systems)
    for k, S in SYS.items():
        log(k, len(S["df"]), f"{S['T'].sum() / 2 / 1e6:.0f} M trips")
    rows = []
    for k, S in SYS.items():
        df, T, Dm = S["df"], S["T"], S["D"]
        ct = np.asarray(df.country.values, dtype=object)
        iu, ju = np.triu_indices(len(df), 1)
        dom = ct[iu] == ct[ju]
        for c in sorted(set(ct)):
            m = dom & (ct[iu] == c)
            if T[iu[m], ju[m]].sum() < 5e6:
                continue
            t, d = T[iu[m], ju[m]], Dm[iu[m], ju[m]]
            rows.append(dict(system=k, market=f"{c} domestic", trips_M=t.sum() / 1e6, gt1500_M=t[d > 1500].sum() / 1e6,
                             gt2500_M=t[d > 2500].sum() / 1e6, gt3500_M=t[d > 3500].sum() / 1e6,
                             paxkm_bn=(t * d).sum() / 1e9))
        t, d = T[iu[~dom], ju[~dom]], Dm[iu[~dom], ju[~dom]]
        if t.sum() > 0:
            rows.append(dict(system=k, market="international (in system)", trips_M=t.sum() / 1e6,
                             gt1500_M=t[d > 1500].sum() / 1e6, gt2500_M=t[d > 2500].sum() / 1e6,
                             gt3500_M=t[d > 3500].sum() / 1e6, paxkm_bn=(t * d).sum() / 1e9))
    pd.DataFrame(rows).to_csv(HERE / "demand_bands.csv", index=False, float_format="%.4g")
    if stage_max < 2:
        return
    BL = {}
    for k, S in SYS.items():
        for dk in ("M2.2", "M2.0"):
            BL[(k, dk)] = cached(f"baseline_{k}_{dk}", lambda: baselines(S, DS[dk]))
        log("baseline", k, len(BL[(k, 'M2.2')]["iu"]))
    if stage_max < 3:
        return
    # points dictionary
    P = dict(WAYPOINTS)
    for k, S in SYS.items():
        for c, la, lo in zip(S["df"].code, S["df"].lat, S["df"].lon):
            P.setdefault(c, (float(la), float(lo)))

    # ------------------------------------------------------------- corridors
    def run_corridors():
        rows, keep = [], {}
        rng = np.random.default_rng(7)
        for name, sysk, chains, note in CORRIDORS:
            edges = chains_to_edges(chains)
            tr = Trunk(edges, P)
            S = SYS[sysk]
            ll = tr.ll
            # climate along spine every ~100 km (subsample)
            ch_pts, ch_brg = [], []
            for a, b, chain in tr.edge_pts:
                c = np.array(chain)[::4]
                for i0, i1 in zip(c[:-1], c[1:]):
                    ch_pts.append(ll[i0]); ch_brg.append(bearing(ll[i0][0], ll[i0][1], ll[i1][0], ll[i1][1]))
            ch_pts, ch_brg = np.array(ch_pts), np.array(ch_brg)
            cl = cutoff_stats(ch_pts, ch_brg, rng)
            gaps = []
            for a, b, chain in tr.edge_pts:
                c = np.array(chain)
                mid = (ll[c[:-1]] + ll[c[1:]]) / 2
                rr = restricted(mid[:, 0], mid[:, 1])
                seg = np.array([tr.G[u][v]["length"] for u, v in zip(c[:-1], c[1:])])
                brk = np.flatnonzero(np.diff(rr.astype(int))) + 1
                for r_, s_ in zip(np.split(rr, brk), np.split(seg, brk)):
                    if not r_[0] and s_.sum() >= 50:
                        gaps.append(f"{a}-{b} {s_.sum():.0f} km")
            for dk, D in DS.items():
                B = BL[(sysk, dk)]
                r = evaluate(tr, S, B, D, R=16, loads=True)
                r_clim = evaluate(tr, S, B, D, R=16, V_bl=cl["mean_eff"] * bt.A_MS * 3.6)
                r_5k = evaluate(tr, S, B, D, R=16, rng_km=RANGE_SENS)
                codes = np.asarray(S["df"].code, dtype=object)
                ph = r["pax"] * r["s_bl"] / 60
                top = np.argsort(-ph)[:8]
                top_s = "; ".join(f"{codes[B['iu'][r['sel'][k]]]}-{codes[B['ju'][r['sel'][k]]]} {ph[k] / 1e6:.2f}" for k in top if ph[k] > 0)
                # share of pax-h from domestic pairs
                same = (np.asarray(S["df"].country.values, dtype=object)[B["iu"][r["sel"]]] ==
                        np.asarray(S["df"].country.values, dtype=object)[B["ju"][r["sel"]]])
                dom_share = float(ph[same].sum() / max(ph.sum(), 1))
                # VOT factor: pax-h weighted income ratio ^0.5 [U]
                ctry = np.asarray(S["df"].country.values, dtype=object)
                sel = r["sel"]
                w = r["pax"] * r["s_bl"]
                gi = np.array([GDPPC.get(c, 20) for c in ctry[B["iu"][sel]]])
                gj = np.array([GDPPC.get(c, 20) for c in ctry[B["ju"][sel]]])
                vf = float((w * np.sqrt(np.sqrt(gi * gj) / 86.0)).sum() / max(w.sum(), 1))
                e = econ(r, B, D, REF)
                ec = econ(r_clim, B, D, REF)
                el = econ(r, B, D, REF, vot_factor=vf)
                peak = max(r["load"].values()) if r["load"] else 0.0
                peak_sub = max(r["load_sub"].values()) if r["load_sub"] else 0.0
                ub_bl = float((B["pax"] * (B["d"] <= RANGE_KM) * np.where(B["t_bl"] - B["t_own"] >= CUT_MIN, B["t_bl"] - B["t_own"], 0)).sum() / 60)
                row = dict(corridor=name, system=sysk, design=dk, km=tr.km, km_overland=tr.km_res,
                           km_overwater=tr.km - tr.km_res, gaps="; ".join(gaps), users_M=r["users_bl"] / 1e6,
                           users_sub_M=r["users_sub"] / 1e6, peak_pax_M=peak / 1e6, peak_pax_sub_M=peak_sub / 1e6,
                           peak_flights_day=peak / (SEATS * LF * 365), paxh_bl_M=r["B_bl"] / 1e6, paxh_sub_M=r["B_sub"] / 1e6,
                           paxh_per_km_bl=r["B_bl"] / max(tr.km_res, 1), system_addressable_paxh_bl_M=ub_bl / 1e6,
                           **e, vot_factor=vf, bc_bl_int_localvot=el["bc_bl_int"],
                           cutoff_share_ge115=cl["share115"], boomless_mach_clim=cl["mean_eff"], Mc_median=cl["Mc_median"],
                           paxh_bl_clim_M=r_clim["B_bl"] / 1e6, bc_bl_int_clim=ec["bc_bl_int"],
                           bc_bl_int_range5000=econ(r_5k, B, D, REF)["bc_bl_int"], domestic_paxh_share=dom_share,
                           top_pairs_Mpaxh=top_s,
                           jet_winter_ms=cl["jet_winter_ms"], jet_u2=cl["jet_u2"], elev_max_km=cl["elev_max"],
                           plateaus="; ".join(plateau_names(ch_pts[:, 0], ch_pts[:, 1])), note=note)
                rows.append(row)
                keep[(name, dk)] = dict(ll=tr.ll, e_u=tr.e_u, e_v=tr.e_v, e_res=tr.e_res, load=r["load"])
                log(f"  {name:45s} {dk} km {tr.km:5.0f} ({tr.km_res:5.0f} overland) users {r['users_bl'] / 1e6:5.1f} M "
                    f"peak {peak / 1e6:5.1f} M paxh {r['B_bl'] / 1e6:5.1f} M  B/C int {e['bc_bl_int']:.3f} net {e['bc_bl_int_net']:.3f} "
                    f"sub {e['bc_sub_int']:.2f} cutoff {cl['share115']:.2f}")
        return pd.DataFrame(rows), keep

    CR, KEEP = cached("corridors_v3", run_corridors)
    CR.to_csv(HERE / "corridors.csv", index=False, float_format="%.4g")
    if stage_max < 4:
        return

    # ------------------------------------------------------------- greedy per scope (M2.2)
    G_rows = []
    for scope, (sysk, groups) in SCOPES.items():
        S = SYS[sysk]
        df = S["df"]
        if groups is None:
            mem = np.ones(len(df), bool)
        else:
            mem = df.country.isin(groups).values | df.bloc.isin(groups).values
        if scope == "US":
            mem = df.code.isin(list(pd.Series(S["T"].sum(1), index=df.code).sort_values(ascending=False).index[:45])).values
        res = cached(f"greedy_{scope}", lambda: greedy(S, BL[(sysk, "M2.2")], DS["M2.2"], P, mem, logf=log))
        if res is None:
            continue
        val, edges, hist = res
        tr = Trunk(edges, P)
        D = DS["M2.2"]
        r = evaluate(tr, S, BL[(sysk, "M2.2")], D, R=16, loads=True)
        e = econ(r, BL[(sysk, "M2.2")], D, REF)
        G_rows.append(dict(scope=scope, edges=" ".join(f"{a}-{b}" for a, b in edges), km=tr.km, km_overland=tr.km_res,
                           users_M=r["users_bl"] / 1e6, peak_pax_M=max(r["load"].values()) / 1e6 if r["load"] else 0,
                           paxh_bl_M=r["B_bl"] / 1e6, paxh_sub_M=r["B_sub"] / 1e6, **e))
        KEEP[(f"greedy {scope}", "M2.2")] = dict(ll=tr.ll, e_u=tr.e_u, e_v=tr.e_v, e_res=tr.e_res, load=r["load"])
        log(f"greedy {scope}: {G_rows[-1]['edges']}  km_res {tr.km_res:.0f}  B/C int {e['bc_bl_int']:.3f}")
        pd.DataFrame(G_rows).to_csv(HERE / "greedy.csv", index=False, float_format="%.4g")
    (CACHE / "keep.pkl").write_bytes(pickle.dumps(KEEP))


# =====================================================================================
# 8. Map
# =====================================================================================
MAP_SET = [  # one spine per family (best or user-requested variant), short label, label offset (dlon, dlat)
    ("US north spine (reference)", "US north (ref)", (-8, 6)), ("US southern spine", "US south", (-6, -7)),
    ("Canada YVR-YYC-YWG-YYZ-YUL", "Canada", (-10, 6)), ("Australia PER-SYD-BNE", "Australia", (-5, 6)),
    ("LON-IST-DXB-BOM-HKG (user)", "LON-IST-DXB-BOM-HKG (user)", (-30, -8)),
    ("China X: PEK-CAN + SHA-CTU", "China X", (6, 4)), ("China NW: PEK-LHW-URC + SHA-XIY-LHW", "China NW", (-12, 7)),
    ("India DEL-HYD-BLR/MAA + CCU", "India", (-24, -14)), ("Russia MOW-SVX-OVB-IKT-KHV-VVO", "Russia", (0, 5)),
    ("Gulf-India-China DXB-DEL-KMG-CAN", "Gulf-India-China", (-2, -24)),
    ("China-SE Asia PEK-CAN-HAN-BKK-SIN", "China-SE Asia", (12, -12)), ("Africa CAI-KRT-ADD-NBO-JNB", "Africa", (5, 0)),
    ("Brazil SAO-BSB + NE + Amazon", "Brazil", (5, 2)), ("S. Cone SCL-BUE-SAO (Andes)", "S. Cone", (-14, -6)),
    ("Andes BOG-LIM-SCL (inland)", "Andes", (-14, 0)), ("Mexico TIJ-GDL-MEX-CUN", "Mexico", (-14, -6)),
]
BINS = [(0.0, 0.05, "#b7d3f6", "< 0.05"), (0.05, 0.2, "#6da7ec", "0.05-0.2"), (0.2, 0.5, "#2a78d6", "0.2-0.5"),
        (0.5, 1.0, "#184f95", "0.5-1.0"), (1.0, 99, "#0d366b", "\u2265 1.0")]


def make_map():
    import json
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    KEEP = pickle.loads((CACHE / "keep.pkl").read_bytes())
    CR = pd.read_csv(HERE / "corridors.csv")
    CR = CR[CR.design == "M2.2"].set_index("corridor")
    gj = json.loads((NM / "routes" / ".cache" / "ne_50m_land.geojson").read_text())
    fig, ax = plt.subplots(figsize=(15, 7.6), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#ffffff")
    for ft in gj["features"]:
        g = ft["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        for pg in polys:
            xy = np.array(pg[0])
            ax.fill(xy[:, 0], xy[:, 1], color="#ecebe7", lw=0, zorder=1)
    for name, lab, off in MAP_SET:
        k = KEEP[(name, "M2.2")]
        bc = float(CR.loc[name, "bc_bl_int"])
        col = [b[2] for b in BINS if b[0] <= bc < b[1]][0]
        ll, load = k["ll"], k["load"]
        for u, v, rr in zip(k["e_u"], k["e_v"], k["e_res"]):
            p = load.get((min(u, v), max(u, v)), 0.0) / 1e6
            x = [ll[u, 1], ll[v, 1]]
            if abs(x[0] - x[1]) > 180:
                continue
            w = 0.6 + 0.11 * p
            ax.plot(x, [ll[u, 0], ll[v, 0]], color=col if rr else "#9a9893", lw=w if rr else max(w * 0.5, 0.6),
                    ls="-" if rr else (0, (2, 1.5)), solid_capstyle="round", zorder=3)
        # label at the most-loaded edge
        e = max(load, key=load.get) if load else (k["e_u"][0], k["e_v"][0])
        x0, y0 = (ll[e[0], 1] + ll[e[1], 1]) / 2, (ll[e[0], 0] + ll[e[1], 0]) / 2
        ax.annotate(f"{lab}\n{CR.loc[name, 'peak_pax_M']:.0f} M pax/yr \u00b7 B/C {bc:.2f}", (x0, y0),
                    (x0 + off[0], y0 + off[1]), fontsize=7.2, color="#2b2b29", zorder=5,
                    arrowprops=dict(arrowstyle="-", color="#9a9893", lw=0.5))
    ax.set_xlim(-130, 160)
    ax.set_ylim(-45, 66)
    ax.set_aspect(1.15)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    h1 = [Line2D([0], [0], color=c, lw=4, label=f"B/C {t}") for _, _, c, t in BINS]
    h2 = [Line2D([0], [0], color="#5f5e5a", lw=0.6 + 0.11 * v, label=f"{v} M pax/yr") for v in (10, 40, 80)]
    h3 = [Line2D([0], [0], color="#9a9893", lw=1.2, ls=(0, (2, 1.5)), label="over water: no shield (design Mach)")]
    ax.legend(handles=h1 + h2 + h3, loc="lower left", fontsize=7.5, frameon=False, ncol=1)
    ax.set_title("BoomShield candidate spines: width = passengers through each km (100 % capture), colour = B/C vs "
                 "Boomless\nM2.2 aircraft, integrated shield cost, gross of fuel/CO\u2082; demand: US = DOT, elsewhere "
                 "calibrated gravity model (unverified)", fontsize=10, color="#2b2b29", loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "global_map.png", facecolor="white")
    log("map written")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=int, default=9)
    ap.add_argument("--figure", action="store_true")
    a = ap.parse_args()
    if a.figure:
        make_map()
    else:
        main(a.stage)
        make_map()
