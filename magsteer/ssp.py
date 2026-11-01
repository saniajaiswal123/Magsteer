"""Sound-speed profiles for the MagSteer simulator. Nootka fault site, NE Pacific."""
import numpy as np


def mackenzie(T, S, z):
    """Mackenzie (1981) nine-term sound speed [m/s]. T degC, S PSU, z m."""
    T = np.asarray(T, float); S = np.asarray(S, float); z = np.asarray(z, float)
    return (1448.96 + 4.591*T - 5.304e-2*T**2 + 2.374e-4*T**3
            + 1.340*(S - 35.0) + 1.630e-2*z + 1.675e-7*z**2
            - 1.025e-2*T*(S - 35.0) - 7.139e-13*T*z**3)


class SSP:
    """Piecewise-linear sound-speed profile."""

    def __init__(self, z, c, name="unnamed", tag="[A]"):
        self.z = np.asarray(z, float)
        self.c = np.asarray(c, float)
        self.name = name
        self.tag = tag
        if np.any(np.diff(self.z) <= 0):
            raise ValueError("depths must strictly increase")

    @property
    def zmax(self):
        return float(self.z[-1])

    def speed(self, depth):
        return np.interp(np.asarray(depth, float), self.z, self.c)

    def axis_depth(self, dz=1.0):
        """Depth of the sound-speed minimum = SOFAR / sound-channel axis."""
        zz = np.arange(0.0, self.zmax + dz, dz)
        return float(zz[int(np.argmin(self.speed(zz)))])

    def mixed_layer_depth(self):
        """Base of the surface duct: first depth where the gradient turns negative."""
        g = np.diff(self.c) / np.diff(self.z)
        neg = np.where(g < 0)[0]
        return float(self.z[neg[0]]) if len(neg) else 0.0

    def resample(self, dz=25.0):
        zz = np.arange(0.0, self.zmax, dz)
        if zz[-1] < self.zmax - 1e-6:
            zz = np.append(zz, self.zmax)   # must land EXACTLY on the bottom depth
        return SSP(zz, self.speed(zz), self.name, self.tag)

    def __repr__(self):
        return (f"SSP({self.name}, 0-{self.zmax:.0f} m, "
                f"axis {self.axis_depth():.0f} m, c_min {self.speed(self.axis_depth()):.1f} m/s)")


# Nootka fault site, ~80 km W of Vancouver Island. Water depth 2372 m [V, Frye 2005].
#
# T/S from WOA23 one-degree ANNUAL objectively-analysed climatology, grid cell
# 49.5 N 128.5 W  [V, NOAA NCEI World Ocean Atlas 2023, woa23_decav_t00mn01 /
# s00mn01]. This REPLACES the previously assumed [A] profile and closes gap #1.
#
# WOA23 terminates at 2200 m at this cell; 2200-2372 m is linearly extrapolated
# in T and S [A, 172 m of a 2372 m column].
_NOOTKA_WOA = [
    (0,    11.021, 32.015), (25,   10.299, 32.230), (50,    8.668, 32.512),
    (75,    7.887, 32.879), (100,   7.496, 33.288), (150,   7.128, 33.769),
    (200,   6.682, 33.897), (300,   5.761, 33.958), (400,   5.143, 34.020),
    (500,   4.728, 34.094), (600,   4.393, 34.165), (700,   4.113, 34.235),
    (800,   3.869, 34.294), (1000,  3.393, 34.384), (1200,  2.975, 34.445),
    (1400,  2.601, 34.493), (1600,  2.311, 34.528), (1800,  2.076, 34.563),
    (2000,  1.902, 34.597), (2200,  1.825, 34.609), (2372,  1.780, 34.620),
]


# Seasonal WOA23 statistical means, same grid cell [V]:
#   woa23_decav_t13/s13 (WINTER), t15/s15 (SUMMER)
# Winter terminates at 2000 m, summer at 2200 m; the final row to 2372 m is
# linearly extrapolated in T and S [A, <16% of the column].
_NOOTKA_SEASONAL = {
    "winter": [
       (0, 8.241, 32.199), (25, 8.191, 32.351), (50, 8.184, 32.475),
       (75, 8.050, 32.914), (100, 7.690, 33.367), (150, 7.311, 33.789),
       (200, 6.881, 33.901), (300, 5.956, 33.958), (400, 5.288, 34.017),
       (500, 4.798, 34.087), (600, 4.461, 34.161), (700, 4.181, 34.235),
       (800, 3.888, 34.293), (1000, 3.392, 34.382), (1200, 2.978, 34.442),
       (1400, 2.605, 34.491), (1600, 2.311, 34.527), (1800, 2.080, 34.560),
       (2000, 1.927, 34.591), (2372, 1.670, 34.651),
    ],
    "summer": [
       (0, 14.467, 31.920), (25, 12.671, 32.187), (50, 8.902, 32.579),
       (75, 7.747, 32.870), (100, 7.369, 33.221), (150, 7.007, 33.754),
       (200, 6.574, 33.898), (300, 5.698, 33.958), (400, 5.109, 34.021),
       (500, 4.714, 34.091), (600, 4.396, 34.159), (700, 4.109, 34.228),
       (800, 3.861, 34.288), (1000, 3.410, 34.380), (1200, 2.989, 34.443),
       (1400, 2.613, 34.491), (1600, 2.317, 34.528), (1800, 2.082, 34.562),
       (2000, 1.904, 34.595), (2372, 1.722, 34.628),
    ],
}


def nootka(season="annual"):
    """
    Nootka fault SSP from WOA23 [V]. season = "annual" | "winter" | "summer".
    All three are NOAA NCEI World Ocean Atlas 2023 statistical means,
    one-degree grid, cell 49.5 N 128.5 W.
    """
    tbl = _NOOTKA_WOA if season == "annual" else _NOOTKA_SEASONAL[season]
    z = np.array([r[0] for r in tbl], float)
    T = np.array([r[1] for r in tbl], float)
    S = np.array([r[2] for r in tbl], float)
    return SSP(z, mackenzie(T, S, z), f"nootka-woa23-{season}",
               "[V] WOA23 49.5N 128.5W")


def from_cast(depth_m, temp_c, sal_psu, name="cast"):
    """Build a [V] profile from any real CTD / WOA / Argo cast."""
    z = np.asarray(depth_m, float)
    return SSP(z, mackenzie(np.asarray(temp_c, float),
                            np.asarray(sal_psu, float), z),
               name, "[V] measured cast")


def isovelocity(zmax=2372.0, c0=1500.0):
    """Straight rays - exact closed form, used for validation."""
    return SSP([0.0, zmax], [c0, c0], "isovelocity", "[V] analytic")
