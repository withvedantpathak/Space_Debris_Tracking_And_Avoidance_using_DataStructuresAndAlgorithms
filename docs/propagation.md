# Orbit propagation (C3 + SGP4)

Code: `src/propagate.[ch]` (Kepler family), `src/sgp4.[ch]` (near-earth SGP4). CLI: `build/dlscms prop <catalog.csv> <dt_s> [flags] [common]`, flags 1 = J2, 2 = drag, 8 = SGP4.
Evidence: `scripts/crosscheck_sgp4.py` -> `results/propagation_error.csv/png`; `scripts/validate_sgp4.py` -> `results/sgp4_validation.csv`.

## Models
| Flags | Model | Used for |
|---|---|---|
| 0 | Kepler two-body | synthetic data (exact ground truth for planted pairs) |
| 1 | + secular J2 (RAAN, argp, mean anomaly; Kozai mean motion handled as in SGP4) | intermediate |
| 3 | + drag proxy: quadratic mean-anomaly term from MEAN_MOTION_DOT | intermediate |
| 8 | **SGP4, hand-written near-earth** (Spacetrack Report #3 / Vallado 2006, WGS-72) | real data |

All are inertial (TEME-like) frames with Earth rotation not applied: screening needs relative geometry only. Epoch is a Julian date in a double (~40 us, ~0.3 m): irrelevant at screening scale.

## SGP4 port: validation
`results/sgp4_validation.csv`, all 19,440 real objects, tsince = 0, 1 h, 1 d, 3 d, 7 d, against the Python `sgp4` package:
- Position difference: max 0.86 mm, median 0.5 mm (the residual is the 1e-6 km print precision of the CLI output). Velocity difference max 0.9 um/s.
- The C code refuses exactly the 805 deep-space objects (period >= 225 min) at 0, 1 h and 1 d. At 3 d and 7 d it refuses a few more (decayed); every one of those is also an error in Python. It never refuses an object Python propagates.
- `tests/test_sgp4.c`: 306 checks over 9 element sets chosen to hit every branch (simplified-drag low perigee, perigee < 98 km and < 156 km, eccentric, equatorial, inclination 180 deg, zero BSTAR), plus error codes, the 225-minute boundary, decay, and the Prop integration.
- Limits: no SDP4; velocity is SGP4's analytic output and differs from d(position)/dt by up to ~1.5 m/s (also true of the reference).

## Kepler-family error against SGP4 (LEO, median km; real catalogue)
| hours after epoch | two-body | + J2 | + J2 + drag | + J2 + drag, 95th pct |
|---|---|---|---|---|
| 0 | 12.0 | 12.0 | 12.0 | 16.1 |
| 1 | 14.6 | 9.9 | 9.9 | 16.6 |
| 6 | 83.7 | 8.4 | 8.5 | 16.2 |
| 24 | 321 | 11.8 | 8.4 | 15.9 |
| 72 | 1235 | 20.3 | 12.4 | 18.2 |
| 168 | 2864 | 85.2 | 11.5 | 27.0 |
| 336 | 5264 | 318 | 17.6 | 189 |

Reading it:
- The ~12 km floor at 0 h is the difference between mean elements and SGP4's osculating state (short-period terms). This is exactly what the SGP4 port removes: with `PROP_SGP4` the model has no such offset.
- J2 is essential (LEO 321 km -> 12 km at 24 h). Drag (confirmed by the BSTAR = 0 experiment) fixes the 95th percentile tail from ~12 h on.
- Large maxima (thousands of km) are decaying objects, where a quadratic drag term is invalid.
- Deep space: 805 objects (GEO/MEO/HEO) have 25 km median error at 0 h and maxima up to 9,000 km. They are outside the project scope (`docs/scope.md`).

## Consequences
- Real-data screening uses SGP4 + LEO scope (flags 12 in `screen`). Event-level agreement with the independent Python implementation: `results/event_agreementscope_sgp4.csv` (12,501 of 12,501 reference events found at 10 km; miss-distance difference median 0.00 km).
- The earlier recommendation to pad thresholds by 25 km (needed for the Kepler model: only 55% of SGP4's 10 km events found at 10 km, 99.7% at 25 km, `results/event_agreement.csv`) no longer applies to the model mismatch. It still applies to the spatial index reach (threshold + v*step) and to altitude-filter bounds.
- Synthetic planted-pair truth is exact only for the two-body model (flags 0).
