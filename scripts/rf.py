"""RF UPLINK: what antenna gain does the satellite link afford?
Produces the -23.2 dBi requirement and the ME antenna comparison."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erfc

# --- Iridium-class L-band link [A: standard link-budget values] ---
F_GHZ, SLANT_KM, P_TX_W, G_SAT_DBK, RATE = 1.616, 1400.0, 7.0, -10.0, 2400.0
EBN0_REQ = 10.5                      # QPSK, BER 1e-6
G_MONOPOLE = 2.0                     # quarter-wave monopole, dBi

FSPL = 20*np.log10(SLANT_KM) + 20*np.log10(F_GHZ) + 92.45
def ebn0(g_dbi):
    eirp = 10*np.log10(P_TX_W) + g_dbi
    return eirp - FSPL + G_SAT_DBK + 228.6 - 10*np.log10(RATE)

print("="*70); print("RF UPLINK BUDGET  (1.616 GHz, 1400 km, 7 W, 2400 bps, QPSK)")
print("="*70)
print("  free-space path loss        %7.2f dB" % FSPL)
print("  Eb/N0 with a %.1f dBi monopole  %7.2f dB" % (G_MONOPOLE, ebn0(G_MONOPOLE)))
print("  QPSK requirement, BER 1e-6  %7.2f dB" % EBN0_REQ)
margin = ebn0(G_MONOPOLE) - EBN0_REQ
g_min  = G_MONOPOLE - margin
print("  MARGIN AVAILABLE            %7.2f dB" % margin)
print("\n  => minimum acceptable antenna gain : %.1f dBi" % g_min)

print("\n" + "="*70); print("PUBLISHED ME ANTENNAS vs THE REQUIREMENT"); print("="*70)
DEV = [("SAW MEMS ME, 430 MHz  [V, Zhang IEEE EDL 2024]", -28.0),
       ("required floor (this work)",                      g_min),
       ("quarter-wave monopole, dry",                      G_MONOPOLE)]
print("  %-46s %8s %10s %9s" % ("device","gain","Eb/N0","verdict"))
print("  " + "-"*78)
for n, g in DEV:
    e = ebn0(g)
    v = "CLOSES" if e >= EBN0_REQ else "%.1f dB short" % (EBN0_REQ - e)
    print("  %-46s %6.1f dBi %7.2f dB %12s" % (n, g, e, v))

print("\n" + "="*70)
print("WHY ME: RESONANCE IS MECHANICAL, NOT ELECTROMAGNETIC")
print("="*70)
print("  ME antenna   f0 = (1/2pi) sqrt(k0/m0)     [V, Wang IET MAP 2024 Eq.1]")
print("               k0 = 14644 N/m, m0 = 0.010 kg -> %.0f Hz (measured 178 Hz) [V]"
      % (np.sqrt(14644/0.010)/(2*np.pi)))
print("               no permittivity term -> no detuning when wetted\n")
print("  conventional f ~ c / (2 L sqrt(eps_r))")
eps_air, eps_sea = 1.0, 80.0
print("               dry  eps_r = %.0f" % eps_air)
print("               wet  eps_r = %.0f  (seawater)" % eps_sea)
print("               resonance shifts by sqrt(80) = %.1fx" % np.sqrt(eps_sea))
print("\n  site sea state [V, DFO C46206, 203,885 records 1988-2022]:")
print("    median Hs 1.95 m | 90th 3.73 m | 95th 4.39 m | 99th 5.73 m")
print("    -> a 0.20 m capsule is intermittently submerged while transmitting")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
g = np.linspace(-40, 5, 400)
ax[0].plot(g, ebn0(g), color="#134074", lw=2.5)
ax[0].axhline(EBN0_REQ, color="#C1443C", ls="--", lw=2)
ax[0].text(-39, EBN0_REQ+1.2, "QPSK requirement, 10.5 dB", fontsize=9, color="#C1443C")
ax[0].axvline(g_min, color="#E8A33D", ls=":", lw=2)
ax[0].text(g_min+0.7, 33, "required floor\n%.1f dBi" % g_min, fontsize=9, color="#B07A20")
ax[0].plot(-28.0, ebn0(-28.0), "o", ms=11, color="#C1443C", mec="k", zorder=5)
ax[0].text(-27, ebn0(-28.0)-4.5, "SAW ME device\n-28 dBi [V]", fontsize=9)
ax[0].plot(G_MONOPOLE, ebn0(G_MONOPOLE), "s", ms=10, color="#2E7D5B", mec="k", zorder=5)
ax[0].fill_between(g, EBN0_REQ, ebn0(g), where=(ebn0(g) >= EBN0_REQ),
                   color="#2E7D5B", alpha=.13)
ax[0].set_xlabel("antenna gain [dBi]"); ax[0].set_ylabel("Eb/N0 [dB]")
ax[0].set_title("Satellite uplink margin\n25.2 dB permits a -23.2 dBi antenna")
ax[0].grid(alpha=.3)

L = np.linspace(0.5, 1.5, 200)
ax[1].plot(L, 1/L, color="#134074", lw=2.5, label="conventional, dry (eps_r = 1)")
ax[1].plot(L, 1/(L*np.sqrt(eps_sea)), color="#C1443C", lw=2.5,
           label="conventional, wetted (eps_r = 80)")
ax[1].axhline(1.0, color="#E8A33D", lw=3, label="ME antenna - mechanical resonance")
ax[1].set_xlabel("normalised electrical length")
ax[1].set_ylabel("resonant frequency, normalised")
ax[1].set_title("Wetting shifts a conventional resonance %.1fx\nthe ME resonance does not move"
                % np.sqrt(eps_sea))
ax[1].set_ylim(0, 2.2); ax[1].grid(alpha=.3); ax[1].legend(fontsize=8.5, loc="upper right")
plt.tight_layout(); plt.savefig("figures/rf_antenna.png", dpi=160)
print("\n  saved -> figures/rf_antenna.png")
