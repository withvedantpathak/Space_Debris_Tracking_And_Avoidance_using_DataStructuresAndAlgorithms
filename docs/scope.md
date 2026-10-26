# Project scope (dataset constraint)

**Decision (DECIDED by the team, 2026-10-21):** the system screens **low-Earth-orbit objects only**: perigee altitude >= 175 km and apogee altitude <= 2000 km above a spherical Earth (R = 6378.137 km). In the CLI this is flag bit 4 of `screen` (`flags` = 4 + 1 for J2 + 2 for drag, e.g. 7).

**Source.** Gomez, Gratl, Bosing, Izzo, "Deterministic Conjunction Tracking in Long-term Space Debris Simulations", 3rd IAA ICSSA, Madrid, arXiv:2203.06957 (PDF read in full, `docs/papers/`). Sec. 3.1: "only satellites operating between 175 and 2000 km above the ground were investigated"; 16,024 objects; sec. 3.2: 10 s time step, objects removed below 150 km; sec. 2: closest approach by linear interpolation inside a step.

**Why this scope also fits our own evidence.**
- Our propagator is accurate to ~12 km median against SGP4 in LEO, but 25 km median (maxima ~9,000 km) for deep-space/HEO objects (`results/propagation_error.csv`). All LEO-scope objects are SGP4 near-earth (period < 225 min).
- Brute force cost is O(n^2) in object count, so the scope also bounds the benchmark.

**Effect on the real catalogue** (CelesTrak active + three debris sets, 19,440 unique objects):
| Filter | Objects kept |
|---|---|
| mean motion >= 6.4 rev/day (earlier near-earth filter, epochs within 3 days of the median) | 18,208 |
| LEO scope 175-2000 km + same epoch and mean-motion filter | 18,191 |
| LEO scope only, no epoch filter | 18,615 (825 excluded: GEO, MEO, HEO, very low perigee) |

The scope removes only 17 more objects than the earlier filter, so the earlier baseline numbers stay valid as a LEO result.

**Deliberate differences from the paper (to state in the report).**
- The paper removes the ISS and the active OneWeb/Starlink satellites (studied separately) and only considers pairs where at least one object is passive or small (< 10 cm). We keep all catalogued objects and all pairs; we have no radius data in the OMM files.
- The paper propagates 20 years with its own integrator; we screen windows of hours to days from mean elements. Our claim is about screening speed and recall, not long-term debris evolution.
- The paper reports no accuracy comparison against SGP4 events; ours is in `results/event_agreement.csv`.

**Out of scope for now.** MEO, GEO, HEO and other deep-space objects (SDP4), objects below 175 km perigee (decaying), manoeuvres, covariance/collision probability.
