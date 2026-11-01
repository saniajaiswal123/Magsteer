"""Bellhop interface: eigenrays -> arrivals. Block 6 of the MagSteer pipeline."""
import os, shutil, subprocess, tempfile
from dataclasses import dataclass
import numpy as np
from .ssp import SSP


def thorpe_alpha(f_khz):
    """Absorption [dB/km]. Khan & Cho Eq. 3 [V]."""
    f2 = f_khz ** 2
    return (0.11*f2/(1+f2) + 44*f2/(4100+f2) + 2.75e-4*f2 + 0.003)


@dataclass
class Arrival:
    amplitude: complex
    delay: float          # s
    angle_src: float      # deg
    angle_rcv: float      # deg  <-- block 7 applies the beam pattern at this angle
    n_surface: int = 0
    n_bottom: int = 0

    @property
    def tl_db(self):
        a = abs(self.amplitude)
        return -20*np.log10(a) if a > 0 else np.inf

    @property
    def path_length(self):
        return self.delay * 1500.0

    def __repr__(self):
        return (f"Arr(TL={self.tl_db:6.2f}dB tau={self.delay*1e3:8.3f}ms "
                f"th_rcv={self.angle_rcv:+7.2f}deg s/b={self.n_surface}/{self.n_bottom})")


class ChannelResponse:
    def __init__(self, arrivals, freq, z_src, z_rcv, range_m):
        self.arrivals = sorted(arrivals, key=lambda a: a.delay)
        self.freq, self.z_src, self.z_rcv, self.range_m = freq, z_src, z_rcv, range_m

    def direct(self):
        d = [a for a in self.arrivals if a.n_surface == 0 and a.n_bottom == 0]
        return max(d, key=lambda a: abs(a.amplitude)) if d else None

    def surface_bounce(self):
        s = [a for a in self.arrivals if a.n_surface == 1 and a.n_bottom == 0]
        return max(s, key=lambda a: abs(a.amplitude)) if s else None

    def is_shadowed(self):
        """No direct arrival = geometry-limited. Nothing for the array to beamform on."""
        return self.direct() is None

    def coherent_tl(self):
        if not self.arrivals: return np.inf
        p = sum(a.amplitude for a in self.arrivals)
        return -20*np.log10(abs(p)) if abs(p) > 0 else np.inf

    def incoherent_tl(self):
        if not self.arrivals: return np.inf
        e = sum(abs(a.amplitude)**2 for a in self.arrivals)
        return -10*np.log10(e) if e > 0 else np.inf

    def delay_spread(self):
        if len(self.arrivals) < 2: return 0.0
        t = np.array([a.delay for a in self.arrivals])
        w = np.array([abs(a.amplitude)**2 for a in self.arrivals]); w /= w.sum()
        return float(np.sqrt((w*(t - (w*t).sum())**2).sum()))

    def __repr__(self):
        return (f"ChannelResponse({len(self.arrivals)} arrivals, "
                f"z={self.z_src:.0f}->{self.z_rcv:.0f} m, r={self.range_m/1000:.2f} km, "
                f"cohTL={self.coherent_tl():.2f} dB)")


def bellhop_path():
    return shutil.which("bellhop.exe") or shutil.which("bellhop")


def eigenrays(ssp: SSP, z_src, z_rcv, range_m, freq,
              n_beams=5000, max_angle=80.0,
              bottom_c=1600.0, bottom_rho=1.8, bottom_alpha=0.5, keep=False):
    """Run Bellhop in arrivals mode and parse the .arr file."""
    exe = bellhop_path()
    if exe is None:
        raise RuntimeError("bellhop.exe not on PATH")

    tmp = tempfile.mkdtemp(prefix="ms_")
    stem = os.path.join(tmp, "case")
    prof = ssp.resample(25.0)
    zb = ssp.zmax

    with open(stem + ".env", "w") as f:
        f.write("'MagSteer'\n")
        f.write(f"{freq:.1f}\n1\n'SVW'\n")
        f.write(f"0 0.0 {zb:.2f}\n")
        for z, c in zip(prof.z, prof.c):
            f.write(f"{z:.2f} {c:.3f} 0.0 1.0 0.0 0.0 /\n")
        f.write("'A' 0.0\n")
        f.write(f"{zb:.2f} {bottom_c:.1f} 0.0 {bottom_rho:.2f} {bottom_alpha:.3f} 0.0 /\n")
        f.write(f"1\n{z_src:.2f} /\n1\n{z_rcv:.2f} /\n1\n{range_m/1000:.6f} /\n")
        f.write(f"'A'\n{n_beams}\n{-max_angle:.2f} {max_angle:.2f} /\n")
        f.write(f"0.0 {zb*1.02:.2f} {range_m/1000*1.02:.6f}\n")

    subprocess.run([exe, stem], cwd=tmp, capture_output=True, timeout=300)

    arr = stem + ".arr"
    if not os.path.exists(arr):
        raise RuntimeError(f"bellhop produced no .arr (see {tmp})")

    rows = []
    with open(arr) as f:
        lines = [l.split() for l in f]
    for p in lines:
        if len(p) < 8:
            continue
        try:
            amp, ph, tau = float(p[0]), float(p[1]), float(p[2])
            a_s, a_r = float(p[4]), float(p[5])
            nt, nb = int(float(p[6])), int(float(p[7]))
        except ValueError:
            continue
        if tau <= 0 or amp <= 0:
            continue
        # Bellhop returns geometric spreading only (alphap = 0 in the .env).
        # Apply Thorpe volume absorption over the actual path length.
        L_km = tau * 1500.0 / 1000.0
        amp *= 10.0 ** (-thorpe_alpha(freq/1000.0) * L_km / 20.0)
        rows.append(Arrival(amp*np.exp(1j*np.radians(ph)), tau, a_s, a_r, nt, nb))

    if not keep:
        shutil.rmtree(tmp, ignore_errors=True)
    return ChannelResponse(rows, freq, z_src, z_rcv, range_m)
