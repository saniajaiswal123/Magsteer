"""FINAL RUN - the three operational metrics. Usage: python3 paper.py [summer|winter]"""
import numpy as np, json, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from magsteer.ssp import nootka
from magsteer.channel import eigenrays
from magsteer.array import LineArray, frye_transducer
from magsteer.link import snr_db, ber_rayleigh, per
from magsteer import constants as C

S = sys.argv[1] if len(sys.argv) > 1 else "summer"
ZB, BW, F = 2372.0, 4000.0, 25000.0
ssp = nootka(S)
RANGES = np.arange(1.0, 8.01, 0.5)
DEPTHS = [400, 600, 800, 1000, 1200, 1400, 1600, 2000]
USABLE = 0.10

def cell(z, r, a):
    try:
        cr = eigenrays(ssp, float(z), ZB, r*1000.0, F, n_beams=8000,
                       bottom_c=C.SEABED_CP_MS,
                       bottom_rho=C.SEABED_DENSITY_KGM3/1000.0,
                       bottom_alpha=C.SEABED_ATTEN_DB_KHZ*F/1e4)
    except Exception:
        return None
    if cr.direct() is None: return None
    g, _ = a.apply(cr); e = 0.0
    for f in np.linspace(F-BW/2, F+BW/2, 25):
        e += abs(sum(x.amplitude*gi*np.exp(-2j*np.pi*(f-F)*x.delay)
                     for x, gi in zip(cr.arrivals, g)))**2
    s = snr_db(C.FRYE_SL_DB, -10*np.log10(e/25), F, a.di_db, BW)
    return float(per(ber_rayleigh(s), C.SEND_PACKET_BYTES))

print(repr(ssp)); print("computing %d ranges..." % len(RANGES))
A = LineArray(n=24, freq=F)
fry, fix, agi, zbest = [], [], [], []
for r in RANGES:
    fry.append(cell(C.FRYE_TRANSDUCER_Z_M, r, frye_transducer(F)))
    fix.append(cell(800.0, r, A))
    c = [(cell(z, r, A), z) for z in DEPTHS]
    c = [(p, z) for p, z in c if p is not None]
    p, z = min(c) if c else (None, None)
    agi.append(p); zbest.append(z)
    print("  %.1f km done" % r, flush=True)

def clean(v): return [1.0 if x is None else x for x in v]
fry, fix, agi = clean(fry), clean(fix), clean(agi)

def radius(p):
    for i in range(len(RANGES)-1):
        if p[i] < USABLE <= p[i+1]:
            f = (USABLE-p[i])/(p[i+1]-p[i])
            return RANGES[i] + f*(RANGES[i+1]-RANGES[i])
    return RANGES[-1] if p[-1] < USABLE else RANGES[0]

print("\n" + "="*72)
print("PACKET ERROR RATE AND RETRANSMISSION COST  (%s, %.0f kHz, SL %.0f dB)"
      % (S, F/1000, C.FRYE_SL_DB))
print("="*72)
print("  %6s | %18s | %18s | %18s | %5s"
      % ("range","Frye 60deg @10m","MagSteer fix 800m","MagSteer agile","z"))
print("  " + "-"*84)
for i, r in enumerate(RANGES):
    tx = lambda p: 1.0/(1.0-p) if p < 0.999 else float('inf')
    print("  %5.1fkm | %8.2f%% %8.2f | %8.2f%% %8.2f | %8.2f%% %8.2f | %4s"
          % (r, fry[i]*100, tx(fry[i]), fix[i]*100, tx(fix[i]),
             agi[i]*100, tx(agi[i]), zbest[i] or "-"))

rf, rx, ra = radius(fry), radius(fix), radius(agi)
print("\n" + "="*72)
print("METRIC 1 - USABLE RADIUS AND AREA SERVED (PER < 10%)")
print("="*72)
for lbl, rr in [("Frye, moored 10 m", rf), ("MagSteer fixed 800 m", rx),
                ("MagSteer depth-agile", ra)]:
    print("  %-24s r = %4.2f km   area = %6.1f km2" % (lbl, rr, np.pi*rr*rr))
print("\n  coverage gain, agile vs Frye : %.2fx area" % ((ra/rf)**2))
print("  relays needed for a 400 km2 array: Frye %.1f  ->  MagSteer %.1f"
      % (400/(np.pi*rf*rf), 400/(np.pi*ra*ra)))

print("\n" + "="*72)
print("METRIC 2 + 3 - ENERGY AND DELIVERY OVER A 13-MONTH DEPLOYMENT")
print("="*72)
print("  Frye: 500 kB/day for 13 months, 1000 bits/joule [V]")
DAYS, KB = 395, 500.0
bits = KB*1024*8*DAYS
print("  payload over the deployment = %.1f MB = %.2e bits" % (KB*DAYS/1024, bits))
print("\n  %6s | %14s | %14s | %14s" % ("range","Frye kJ","agile kJ","energy saved"))
print("  " + "-"*62)
for i, r in enumerate(RANGES):
    if r not in (2.0, 3.0, 4.0, 5.0, 6.0): continue
    ef = bits*(1.0/(1.0-fry[i]))/1000.0/1000.0 if fry[i] < 0.999 else float('inf')
    ea = bits*(1.0/(1.0-agi[i]))/1000.0/1000.0 if agi[i] < 0.999 else float('inf')
    sv = "n/a" if not np.isfinite(ef) else ("%.1f%%" % (100*(ef-ea)/ef))
    print("  %5.1fkm | %13.1f | %13.1f | %13s" % (r, ef, ea, sv))

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].semilogy(RANGES, np.array(fry)*100, "s-", color="#999999", lw=2, ms=6,
               label="Frye, moored 10 m, 60 deg")
ax[0].semilogy(RANGES, np.array(fix)*100, "^-", color="#1f78b4", lw=2, ms=6,
               label="MagSteer, fixed 800 m")
ax[0].semilogy(RANGES, np.array(agi)*100, "o-", color="#d95f02", lw=2.5, ms=7,
               label="MagSteer, depth-agile")
ax[0].axhline(10, color="k", ls="--", lw=1.5)
ax[0].text(1.1, 11.5, "10% PER - usable link", fontsize=9)
ax[0].set_xlabel("OBS node horizontal range [km]")
ax[0].set_ylabel("packet error rate [%]")
ax[0].set_title("Link quality vs range\nNootka fault, %s, %.0f kHz" % (S, F/1000))
ax[0].grid(alpha=.3, which="both"); ax[0].legend(fontsize=9)

tx = lambda v: [1.0/(1.0-p) if p < 0.99 else 100 for p in v]
ax[1].plot(RANGES, tx(fry), "s-", color="#999999", lw=2, ms=6, label="Frye")
ax[1].plot(RANGES, tx(agi), "o-", color="#d95f02", lw=2.5, ms=7, label="MagSteer agile")
ax[1].set_ylim(0.9, 5); ax[1].axhline(1.33, color="g", ls=":", lw=2)
ax[1].text(1.1, 1.38, "Send et al. deployed: 1.33 tx/packet [V]", fontsize=9, color="g")
ax[1].set_xlabel("OBS node horizontal range [km]")
ax[1].set_ylabel("transmissions per delivered packet")
ax[1].set_title("Retransmission cost\n(lower = less energy, longer battery)")
ax[1].grid(alpha=.3); ax[1].legend(fontsize=9)
plt.tight_layout(); plt.savefig("figures/paper_%s.png" % S, dpi=160)
json.dump(dict(ranges=RANGES.tolist(), frye=fry, fixed=fix, agile=agi,
               zbest=[z or 0 for z in zbest]),
          open("out/paper_%s.json" % S, "w"), indent=1)
print("\n  saved -> figures/paper_%s.png  +  out/paper_%s.json" % (S, S))
