"""Link budget and error rates. Khan & Cho equations [V]."""
import numpy as np


def noise_level_db(f_hz, bw_hz=4000.0):
    """NL = 50 - 18 log10(f_kHz) dB re uPa/Hz, integrated over bandwidth [V, Khan & Cho]."""
    return 50.0 - 18.0*np.log10(f_hz/1000.0) + 10.0*np.log10(bw_hz)


def snr_db(sl_db, tl_db, f_hz, di_db=0.0, bw_hz=4000.0):
    """SNR = SL - TL - NL + DI [V, Khan & Cho Eq. 1]."""
    return sl_db - tl_db - noise_level_db(f_hz, bw_hz) + di_db


def ber_rayleigh(snr_db_val):
    """BPSK in Rayleigh fading: Pb = 0.5[1 - sqrt(g/(1+g))] [V, Khan & Cho Eq. 7]."""
    g = 10.0**(np.asarray(snr_db_val, float)/10.0)
    return 0.5*(1.0 - np.sqrt(g/(1.0+g)))


def per(ber, packet_bits=1024):
    """PER = 1 - (1-Pb)^L [V, Khan & Cho Eq. 8]."""
    return 1.0 - (1.0 - np.asarray(ber, float))**packet_bits


def snr_from_per(target_per, packet_bits=1024):
    """Invert PER -> SNR. Used to back out the glider's operating point from its measured 75%."""
    pb = 1.0 - (1.0 - target_per)**(1.0/packet_bits)
    g = (1.0 - 2*pb)**2 / (1.0 - (1.0 - 2*pb)**2)
    return 10.0*np.log10(g)


def shannon_bps(snr_db_val, bw_hz=4000.0):
    return bw_hz*np.log2(1.0 + 10.0**(np.asarray(snr_db_val, float)/10.0))
