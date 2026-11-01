"""RF DATA RATE TRADE: link margin vs surface time."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

F, S, P, GSAT, EBREQ = 1.616, 1400.0, 7.0, -10.0, 10.5
FSPL = 20*np.log10(S) + 20*np.log10(F) + 92.45
def ebn0(g, R): return 10*np.log10(P) + g - FSPL + GSAT + 228.6 - 10*np.log10(R)

RATES = [2400, 1200, 600, 300, 150]
G_ME, G_MONO, KB = -28.0, 2.0, 100.0

print("="*72)
print("RF DATA RATE TRADE  -  ME antenna -28 dBi [V, Zhang IEEE EDL 2024]")
print("="*72)
print("  %7s | %9s | %9s | %10s | %14s"
      % ("rate", "Eb/N0", "margin", "verdict", "surface time"))
print("  " + "-"*64)
for R in RATES:
    e = ebn0(G_ME, R); m = e - EBREQ
    t = KB*1024*8/R/60
    print("  %6dbps | %6.2f dB | %6.2f dB | %10s | %8.1f min"
          % (R, e, m, "CLOSES" if m >= 0 else "%.1f short" % -m, t))

print("\n  payload assumed: %.0f kB per surfacing" % KB)
print("  600 bps is Send et al.'s own usable acoustic baud [V]")
print("\n  DESIGN POINT: 600 bps closes the link with 1.2 dB spare,")
print("                at the cost of %.0f min surfaced instead of %.0f min."
      % (KB*1024*8/600/60, KB*1024*8/2400/60))

fig, ax = plt.subplots(figsize=(7.5, 5))
r = np.logspace(np.log10(100), np.log10(4800), 300)
ax.semilogx(r, ebn0(G_ME, r), color="#134074", lw=2.5,
            label="ME antenna, -28 dBi [V]")
ax.semilogx(r, ebn0(G_MONO, r), color="#8DA9C4", lw=2, ls="--",
            label="quarter-wave monopole, +2 dBi")
ax.axhline(EBREQ, color="#C1443C", ls="--", lw=2)
ax.text(110, EBREQ+0.8, "QPSK requirement, 10.5 dB", fontsize=9, color="#C1443C")
ax.axvline(600, color="#E8A33D", ls=":", lw=2.5)
ax.text(640, 30, "design point\n600 bps", fontsize=10, color="#B07A20", weight="bold")
ax.fill_between(r, EBREQ, ebn0(G_ME, r), where=(ebn0(G_ME, r) >= EBREQ),
                color="#2E7D5B", alpha=.15)
for R in RATES:
    ax.plot(R, ebn0(G_ME, R), "o", ms=8, color="#134074", mec="k", zorder=5)
ax.set_xlabel("uplink data rate [bps]"); ax.set_ylabel("Eb/N0 [dB]")
ax.set_title("A published ME antenna closes the satellite link at 600 bps\n"
             "(shaded = link closes)")
ax.grid(alpha=.3, which="both"); ax.legend(fontsize=9, loc="upper right")
plt.tight_layout(); plt.savefig("figures/rf_rate.png", dpi=160)
print("\n  saved -> figures/rf_rate.png")
