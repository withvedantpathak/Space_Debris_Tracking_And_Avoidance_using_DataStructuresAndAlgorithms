# Synthetic catalogue model (C2)

Generator: `src/synth.c`; CLI: `build/dlscms gen <n> <seed> <out.csv> [pairs]`. Same seed gives the same file on the same machine/compiler (the draws use libm `sin/cos/log/exp`, so the last printed digit could in principle differ on another platform). Output is OMM CSV in the CelesTrak column layout, so it goes through the same parser as real data.

## Filler objects: where the numbers come from
Histograms taken from the CelesTrak snapshot in `data/raw/` (active + COSMOS 2251 + FENGYUN 1C + IRIDIUM 33 debris, 19,440 unique NORAD ids, epochs 2026-09-21..2026-10-22). Counts are hard-coded in `synth.c`.

| Property | Model | Source / caveat |
|---|---|---|
| Orbit class | 95.8% LEO (a-R < 2000 km), 4.2% above (814/19,440) | measured |
| LEO altitude | 100 km bins from 0 km, weights = measured counts, uniform within bin, floor 200 km | measured; dominated by the 400-500 km bin (10,618 of 18,626) |
| Above LEO | 70% GEO (35,786 km +-150), 30% uniform 2,000-25,000 km | coarse split chosen by hand from percentiles; not fitted |
| Inclination | 10 degree bins, measured weights, uniform within bin | measured, LEO objects only; applied to all classes; no altitude-inclination correlation |
| Eccentricity | lognormal, median 2e-4, sigma 2.56, capped at 0.05 and at perigee >= 200 km | matched to LEO median 0.0002 and 90th percentile 0.0053; real 99th percentile is 0.033 |
| RAAN, arg. of perigee, mean anomaly | uniform 0-360 | simplification |
| Epoch | all objects at 2026-10-20T00:00:00 | real epochs differ by days |
| BSTAR, drag terms | 0 | not modelled |

Known mismatch: median LEO inclination is 59 deg in the synthetic set versus 53 deg in the real catalogue (real objects pile up near 53 deg, synthetic ones are spread uniformly inside each 10 deg bin; `results/dataset_compare.png`, `results/dataset_stats.csv`).

Not modelled: HEO/Molniya population (real catalogue has 45 objects with e >= 0.25), orbital shells/constellation structure, clustering of debris from breakups. Real data is used for the final scaling runs.

## Planted pairs (known answers)
Two circular orbits (e = 0) of equal radius (450-1000 km altitude), inclinations 40-140 degrees, planes differing by 15-60 degrees, crossing angle between 20 and 150 degrees. B is shifted along its track until the exact two-body closest approach equals the designed miss distance (golden-section search on the true circular-orbit motion, iterated to 1e-7 km). Closest approach is 0.5-20 hours after epoch. Designed miss distances cycle through 0.5, 1, 2, 4, 6, 8 km (conjunctions at the default 10 km threshold) and 12, 15, 25 km (near misses that must NOT be reported at 10 km).

`<name>_truth.csv`: `id_a,id_b,tca_s_after_epoch,miss_km,rel_speed_kms`. Truth is exact for the two-body model only; it is not claimed for SGP4 or J2. `tests/test_catalog.c` checks every planted pair with a separate propagator written in the test file, using the elements re-read from the CSV.

Other (random) objects can also come close to each other; the truth file lists only planted pairs, so use the brute-force baseline (C4) as the reference for total counts and these pairs as a recall floor.

## Datasets in `data/synthetic/`
`synth_1000`, `synth_5000`, `synth_20000` (seed 1; 10 / 50 / 100 planted pairs). Larger sets: `build/dlscms gen 50000 1 data/synthetic/synth_50000.csv`.
