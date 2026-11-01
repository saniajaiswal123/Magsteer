"""Verified constants. Every value is [V] with its source."""

# Frye et al. 2005, Proc. OCEANS 2005 Europe 969-974, DOI 10.1109/OCEANSE.2005.1513188
FRYE_SL_DB          = 190.0    # "The source level is 190 dB"                [V]
FRYE_TX_POWER_W     = 8.0      # "transmit power is approximately 8 Watts"   [V]
FRYE_BEAMWIDTH_DEG  = 60.0     # "directional transducers with 60 deg beams" [V]
FRYE_TRANSDUCER_Z_M = 10.0     # bottom of a 10 m chain below the u-joint    [V]
FRYE_WATER_DEPTH_M  = 2362.0   # "deployed in 2362 m of water"               [V]
FRYE_RANGE_KM       = 4.0      # "longer range (~4km) acoustic telemetry"    [V]
FRYE_RATE_BPS       = 5400.0   # "higher data rate (~5400bps)"               [V]
FRYE_BITS_PER_JOULE = 1000.0   # "on order 1000 bits per joule"              [V]
FRYE_MAX_SENSORS    = 15       # "Multiple remote sensors (up to 15)"        [V]
FRYE_KB_PER_DAY     = 500.0    # ">500 Kbytes/day from two instruments"      [V]
FRYE_COVERAGE_KM2   = (5.0, 20.0)      # "sensors over a 5-20 km^2 area"     [V]
FRYE_CTRL_BPS       = 80.0             # FH-FSK control                      [V]
FRYE_DATA_BPS_RANGE = (300.0, 5300.0)  # PSK data uplink                     [V]

# Send, Regier & Jones 2013, JAOT 30(5) 984-998, DOI 10.1175/JTECH-D-11-00169.1
SEND_SUCCESS_RATE   = 0.75     # "success rate of 75%"                       [V]
SEND_J_PER_BYTE     = 0.49     # measured field EEPB                         [V]
SEND_J_PER_BYTE_MOD = 0.21     # model prediction                            [V]
SEND_OMNI_SL_DB     = 185.0    # omnidirectional transducer SL               [V]
SEND_DIR_SL_DB      = 190.0    # directional transducer SL                   [V]
SEND_RANGE_KM       = (3.5, 5.5)          # typical horizontal range         [V]
SEND_RANGE_MAX_KM   = 7.0                 # max, 8 km slant                  [V]
SEND_MODEM_BAND_HZ  = (9000.0, 14000.0)   # Benthos ATM-885 LF               [V]
SEND_USABLE_BAUD    = 600.0    # "only 600 and lower were usable"            [V]
SEND_PACKET_BYTES   = 1024     # errors rise sharply beyond 1024 bytes       [V]
# Send Fig. 1: source at 1000 m -> shadow zone beyond 6 km;
#              source at 3200 m -> no near-surface shadow zone.              [V]

# Hamilton 1980, JASA 68(5) 1313-1340, via Dong et al. 2012 OJA 2, 25-33
SEABED_DENSITY_KGM3 = 1767.0   # silt                                        [V]
SEABED_CP_MS        = 1623.0   # compressional velocity                      [V]
SEABED_ATTEN_DB_KHZ = 0.673    # compressional attenuation                   [V]

# THK ball screw product information
BALLSCREW_ETA       = 0.90     # positive efficiency                         [V]
SLIDING_SCREW_ETA   = 0.32     # sliding screw, for comparison               [V]
