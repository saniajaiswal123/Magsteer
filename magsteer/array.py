"""Vertical line array: steering, beam pattern, per-arrival gain. Block 7."""
import numpy as np

C0 = 1500.0


class LineArray:
    """
    Uniform vertical line array, N elements at lambda/2.

    Vertical because the direct path and its surface reflection are separated
    in ELEVATION - so elevation resolution is what buys multipath rejection.
    """

    def __init__(self, n=24, freq=12000.0, spacing_frac=0.5, c=C0):
        self.n, self.freq, self.c = n, freq, c
        self.lam = c / freq
        self.d = spacing_frac * self.lam
        self.aperture = (n - 1) * self.d

    @property
    def di_db(self):
        """Directivity index, 10 log10(N). Validated by Zhao: N=2 -> 2.69 dB measured vs 3.01 predicted [V]."""
        return 10.0 * np.log10(self.n)

    @property
    def beamwidth_deg(self):
        """3 dB beamwidth at broadside, uniform weighting."""
        return np.degrees(0.886 * self.lam / (self.n * self.d))

    def steering_vector(self, theta_deg):
        k = 2 * np.pi * self.freq / self.c
        m = np.arange(self.n)
        return np.exp(1j * k * m * self.d * np.sin(np.radians(theta_deg)))

    def response(self, theta_deg, steer_deg=0.0):
        """Normalised complex array response (1.0 on boresight)."""
        th = np.atleast_1d(np.asarray(theta_deg, float))
        w = self.steering_vector(steer_deg)
        a = np.exp(1j * 2*np.pi*self.freq/self.c *
                   np.arange(self.n)[None, :] * self.d * np.sin(np.radians(th))[:, None])
        r = (a @ np.conj(w)) / self.n
        return r[0] if np.isscalar(theta_deg) or np.ndim(theta_deg) == 0 else r

    def gain_db(self, theta_deg, steer_deg=0.0):
        r = np.abs(self.response(theta_deg, steer_deg))
        return 20 * np.log10(np.maximum(r, 1e-12))

    def apply(self, channel, steer_deg=None):
        """
        Per-arrival complex gains. THIS is where beamforming becomes physical:
        the direct path sits in the mainlobe, the surface bounce lands wherever
        geometry puts it and takes the sidelobe level at that angle.
        """
        if steer_deg is None:
            d = channel.direct()
            steer_deg = d.angle_rcv if d else 0.0
        return np.array([self.response(a.angle_rcv, steer_deg) for a in channel.arrivals]), steer_deg

    def __repr__(self):
        return (f"LineArray(N={self.n}, f={self.freq/1000:.0f} kHz, d={self.d*100:.2f} cm, "
                f"aperture={self.aperture:.2f} m, DI={self.di_db:.2f} dB, "
                f"beamwidth={self.beamwidth_deg:.2f} deg)")


class ConicalBeam:
    """
    Fixed conical-beam transducer, as used by Frye et al. 2005.

    "The modems all use directional transducers with 60 degree beams... The
    source level is 190 dB and the transmit power is approximately 8 Watts."
    [V, Frye et al. 2005, sec. II-F]
    """

    def __init__(self, beamwidth_deg=60.0, freq=25000.0):
        self.beamwidth = beamwidth_deg
        self.freq = freq
        half = np.radians(beamwidth_deg / 2.0)
        self.omega = 2 * np.pi * (1 - np.cos(half))

    @property
    def n(self):
        return 1

    @property
    def di_db(self):
        return 10.0 * np.log10(4 * np.pi / self.omega)

    @property
    def beamwidth_deg(self):
        return self.beamwidth

    def gain_db(self, theta_deg, steer_deg=0.0):
        """Flat inside the cone, -20 dB outside: a fixed 60 deg beam does not
        discriminate within its own width."""
        d = np.abs(np.asarray(theta_deg, float) - steer_deg)
        return np.where(d <= self.beamwidth / 2.0, 0.0, -20.0)

    def response(self, theta_deg, steer_deg=0.0):
        return 10.0 ** (self.gain_db(theta_deg, steer_deg) / 20.0)

    def apply(self, channel, steer_deg=None):
        if steer_deg is None:
            d = channel.direct()
            steer_deg = d.angle_rcv if d else 0.0
        g = np.array([complex(self.response(a.angle_rcv, steer_deg))
                      for a in channel.arrivals])
        return g, steer_deg

    def __repr__(self):
        return f"ConicalBeam({self.beamwidth:.0f} deg, DI={self.di_db:.2f} dB) [V, Frye 2005]"


def frye_transducer(freq=25000.0):
    """Frye et al. 2005 baseline: 60 deg beam, DI 11.74 dB [V]."""
    return ConicalBeam(60.0, freq)


def omni():
    """Send et al. glider: 'nearly omnidirectional' transducer, DI = 0 [V]."""
    return LineArray(n=1, freq=25000.0)
