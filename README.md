# MagSteer

A depth-agile vertical array for bidirectional acoustic data recovery from
ocean-bottom seismometers.

Ocean-bottom seismometers record for a year or more, but their data stays
inaccessible until a vessel recovers the instrument. Acoustically linked
surface moorings solve this within a few kilometres. Their reach is limited
by an effect transmit power cannot fix: a signal from a seafloor node reaches
the relay twice, directly and by surface reflection, and the copies interfere.

This repository contains the model that quantifies that limit and the
architecture that addresses it.

## Headline results

| | Moored baseline | Depth-agile array |
|---|---|---|
| Surface-echo suppression | under 3 dB | 13.6 - 32.9 dB |
| Usable radius, summer | 3.22 km | 5.22 km |
| Usable radius, winter | 4.51 km | 5.21 km |
| Seasonal variation | 40% | 0.2% |
| Energy, 193 MB at 6 km summer | 5,625.7 kJ | 2,113.6 kJ |

An ablation shows the same array held at a fixed depth reaches only 3.07 km
in winter, below the moored baseline. The result follows from combining
angular resolution with depth selection, not from either alone.

## Running it

    python3 -m venv venv && source venv/bin/activate
    pip install -r requirements.txt

Bellhop is not a Python package. Build it from the Acoustics Toolbox
(oalib-acoustics.org) and put it on your PATH:

    export PATH="$HOME/at/bin:$PATH"
    which bellhop.exe

Then:

    python3 scripts/validate.py
    python3 scripts/paper.py summer
    python3 scripts/paper.py winter
    python3 scripts/rf.py

Each script writes JSON to out/ and figures to figures/.

## Layout

    magsteer/constants.py   every published value with its source quoted
    magsteer/ssp.py         WOA23 profiles, Mackenzie equation
    magsteer/channel.py     Bellhop interface, Thorpe absorption per path
    magsteer/array.py       LineArray, ConicalBeam baseline, omni
    magsteer/link.py        Khan and Cho link budget, Rayleigh BER, PER
    scripts/                run scripts producing the published figures
    paper/                  conference and journal manuscripts

## Data provenance

Every environmental input is measured. Nothing was assumed or chosen for
convenience.

| Input | Source |
|---|---|
| Sound speed | NOAA World Ocean Atlas 2023, 1-deg means at 49.5N 128.5W |
| Sound-speed equation | Mackenzie nine-term, JASA 70(3), 1981 |
| Bathymetry | GEBCO 2026 Grid, 2341 m at site, 0.2 m/km slope |
| Wave climate | Fisheries and Oceans Canada C46206, 203,885 records 1988-2022 |
| Seabed | Hamilton, JASA 68(5), 1980: 1767 kg/m3, 1623 m/s, 0.673 dB/kHz |
| Propagation | Bellhop, Acoustics Toolbox (Porter) |
| Baseline hardware | Frye et al., OCEANS 2005: 60 deg beam, 190 dB, 10 m |
| Mobile benchmark | Send, Regier and Jones, JTECH 30(5), 2013 |

WOA23 and GEBCO files are not committed; they are large and freely available
from source.

## Validation

The ray engine reproduces the closed form for straight rays in an isovelocity
ocean to 0.0000 degrees and 0.00 m. More usefully, the model reproduces the
published operating range of two independently deployed systems without being
fitted to either: Frye fails beyond 4 km and reports about 4 km; the Send
glider degrades beyond 3 km and reports 3.5-5.5 km typical.

## Limitations

This work is simulation. It has not been validated against tank or sea data.

- Bellhop is a 2-D ray model over a horizontally uniform ocean: no currents,
  internal waves, roughness, or out-of-plane propagation.
- The array is steered exactly at the direct arrival. A real system estimates
  that angle with error; on a 4.23 degree beam, one or two degrees costs
  several dB of main-lobe gain. Adaptive weighting is not implemented.
- WOA23 ends at 2200 m here; the final 172 m is linearly extrapolated.
- The 60 degree baseline is modelled flat-top, so suppression is reported as
  a bound rather than a value.
- Beyond about 15 km the two arrivals converge below the beamwidth.

## Known issues

constants.py defines SEND_PACKET_BYTES = 1024, passed to per(ber,
packet_bits). The value is correct and traces to Khan and Cho's 1024-bit
packet, but the name and attribution are wrong.

## License

MIT for the code. The manuscripts are not covered by it; rights may transfer
to the publisher on acceptance.
