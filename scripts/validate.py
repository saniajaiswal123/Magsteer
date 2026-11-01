import numpy as np
from magsteer.ssp import isovelocity
from magsteer.channel import eigenrays, thorpe_alpha

zs, zr, r, f = 200.0, 2371.0, 3000.0, 12000.0
L  = np.hypot(r, zr - zs)
th = np.degrees(np.arctan2(zr - zs, r))
Ls = np.hypot(r, zr + zs)
ths = np.degrees(np.arctan2(zr + zs, r))

print("VALIDATION: isovelocity, straight rays, exact closed form")
print(f"  analytic direct : L={L:8.1f} m  theta_rcv={th:+7.3f} deg  tau={L/1500:.6f} s")
print(f"  analytic surface: L={Ls:8.1f} m  theta_rcv={ths:+7.3f} deg  tau={Ls/1500:.6f} s")

cr = eigenrays(isovelocity(), zs, zr, r, f)
d, s = cr.direct(), cr.surface_bounce()
print(f"\n  bellhop: {len(cr.arrivals)} arrivals")
print(f"  direct  : L={d.path_length:8.1f} m  theta_rcv={d.angle_rcv:+7.3f} deg  "
      f"err {d.angle_rcv-th:+.4f} deg / {d.path_length-L:+.2f} m")
if s:
    print(f"  surface : L={s.path_length:8.1f} m  theta_rcv={s.angle_rcv:+7.3f} deg  "
          f"err {abs(s.angle_rcv)-ths:+.4f} deg / {s.path_length-Ls:+.2f} m")
    print(f"  sep     : {abs(s.angle_rcv - d.angle_rcv):.2f} deg   "
          f"amp ratio {abs(s.amplitude)/abs(d.amplitude):.3f}")
print(f"\n  alpha(12 kHz) = {thorpe_alpha(12.0):.3f} dB/km")
print(f"  delay spread  = {cr.delay_spread()*1e3:.3f} ms")
print("\n  strongest arrivals:")
for a in sorted(cr.arrivals, key=lambda x: -abs(x.amplitude))[:5]:
    print("   ", a)
