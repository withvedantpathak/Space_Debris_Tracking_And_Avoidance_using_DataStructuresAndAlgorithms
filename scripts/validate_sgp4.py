"""Validate src/sgp4.c against the Python `sgp4` package (reference implementation, WGS-72) on the whole real catalogue.
For every object, propagate tsince = 0, 1 h, 1 d, 3 d, 7 d from its own element epoch with both and compare TEME position (km) and
velocity (km/s). Objects that the C code refuses (deep space) are reported separately and must be exactly the deep-space ones.
Writes results/sgp4_validation.csv. Exit code 1 if any near-earth object differs by more than 1e-6 km (1 mm, the CSV print precision) at any time or if C refuses an object Python propagates."""
import os, subprocess, io, sys, numpy as np, pandas as pd
from sgp4 import omm
from sgp4.api import Satrec
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
src = os.path.join(root, "build", "merged_real.csv")
exe = os.path.join(root, "build", "dlscms")
with open(src) as fh: rows = list(omm.parse_csv(fh))
sats = {}
for r in rows:
    s = Satrec(); omm.initialize(s, r); sats[int(r["NORAD_CAT_ID"])] = s
out = []; worst = 0.0
for label, dt in [("0", 0), ("1h", 3600), ("1d", 86400), ("3d", 259200), ("7d", 604800)]:
    c = pd.read_csv(io.StringIO(subprocess.run([exe, "prop", src, str(dt), "11"], capture_output=True, text=True, check=True).stdout))
    # flags 11 = J2 + drag + SGP4; objects the C code could not propagate are parked at x >= 1e9
    parked = c.x_km >= 1e8
    dr, dv, nref_err, n = [], [], 0, 0
    for row in c[~parked].itertuples():
        e, r, v = sats[row.norad_id].sgp4_tsince(dt / 60.0)
        if e != 0: nref_err += 1; continue
        dr.append(np.linalg.norm(np.array(r) - [row.x_km, row.y_km, row.z_km])); dv.append(np.linalg.norm(np.array(v) - [row.vx_kms, row.vy_kms, row.vz_kms]))
    dr, dv = np.array(dr), np.array(dv)
    deep = {i for i, s in sats.items() if s.method == 'd'}
    park_ids = set(c.norad_id[parked])
    # every object C refused that is not deep space must also be an error in Python at this time (decayed etc.)
    c_only = sum(1 for i in park_ids - deep if sats[i].sgp4_tsince(dt / 60.0)[0] == 0)
    out.append(dict(t=label, compared=len(dr), c_refused=int(parked.sum()), c_refused_equals_deep_space=(park_ids == deep), c_refused_but_python_ok=c_only,
                    python_error_where_c_ok=nref_err, max_pos_diff_km=dr.max(), p99_pos_diff_km=np.percentile(dr, 99), median_pos_diff_km=np.median(dr),
                    max_vel_diff_kms=dv.max()))
    worst = max(worst, dr.max())
res = pd.DataFrame(out); res.to_csv(os.path.join(root, "results", "sgp4_validation.csv"), index=False)
pd.set_option("display.width", 250); print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
ok = worst <= 1e-6 and (res.c_refused_but_python_ok == 0).all()
print("RESULT:", "PASS" if ok else "FAIL", f"(worst position difference {worst:.3g} km)")
sys.exit(0 if ok else 1)
