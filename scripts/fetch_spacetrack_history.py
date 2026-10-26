"""L34 / F10: download element-set HISTORY from Space-Track (gp_history) for a reproducible sample of the project's LEO objects,
so that scripts/calibrate_age_from_history.py can measure how fast prediction error grows with element age.

STATUS: written but NOT tested against the live service (no account was available when it was written). Run with --dry-run first.
Credentials: ONLY from the environment variables SPACETRACK_USER and SPACETRACK_PASS. They are never written to a file or printed.
Politeness: Space-Track asks for fewer than 30 requests per minute and 300 per hour; this script waits >= 6 s between requests
(max 10 per minute) and puts many NORAD ids into one request. Check Space-Track's terms of use before sharing the downloaded file:
it is written to data/history/ (not part of the delivered data) and should not be redistributed.

Usage: python3 scripts/fetch_spacetrack_history.py [--dry-run] [--objects 300] [--days 7] [--batch 50] [--seed 1]
Output: data/history/spacetrack_gp_history.csv  (OMM CSV, several rows per NORAD_CAT_ID)
Then:   python3 scripts/calibrate_age_from_history.py data/history/spacetrack_gp_history.csv 5 95"""
import argparse, csv, http.cookiejar, os, random, sys, time, urllib.parse, urllib.request

ap = argparse.ArgumentParser()
ap.add_argument("--dry-run", action="store_true", help="print the request URLs and exit; no network, no credentials needed")
ap.add_argument("--objects", type=int, default=300, help="sample size (random, seeded) of LEO-scope NORAD ids")
ap.add_argument("--days", type=int, default=7, help="history window: elements with EPOCH in the last N days")
ap.add_argument("--batch", type=int, default=50, help="NORAD ids per request")
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
BASE = "https://www.space-track.org"
MU, RE = 398600.4418, 6378.137

ids = {}
for fn in ("active", "cosmos-2251-debris", "fengyun-1c-debris", "iridium-33-debris"):
    with open(os.path.join(ROOT, "data", "raw", fn + ".csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            n = float(r["MEAN_MOTION"]) * 2 * 3.141592653589793 / 86400; ecc = float(r["ECCENTRICITY"])
            sma = (MU / n ** 2) ** (1 / 3)
            if sma * (1 - ecc) - RE >= 175 and sma * (1 + ecc) - RE <= 2000: ids[int(r["NORAD_CAT_ID"])] = 1   # project scope
pool = sorted(ids); random.Random(a.seed).shuffle(pool); pick = sorted(pool[:a.objects])
batches = [pick[i:i + a.batch] for i in range(0, len(pick), a.batch)]
def url(b):
    q = f"/basicspacedata/query/class/gp_history/NORAD_CAT_ID/{','.join(map(str, b))}/EPOCH/>now-{a.days}/orderby/NORAD_CAT_ID asc,EPOCH asc/format/csv"
    return BASE + urllib.parse.quote(q, safe="/,>:")
print(f"{len(pick)} objects in {len(batches)} requests, window {a.days} days")
if a.dry_run:
    for b in batches: print(url(b))
    sys.exit(0)

user, pw = os.environ.get("SPACETRACK_USER"), os.environ.get("SPACETRACK_PASS")
if not user or not pw: sys.exit("set SPACETRACK_USER and SPACETRACK_PASS in the environment first (do not put them in a file)")
op = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
resp = op.open(urllib.request.Request(BASE + "/ajaxauth/login", data=urllib.parse.urlencode({"identity": user, "password": pw}).encode()), timeout=60)
if b"Failed" in resp.read(): sys.exit("login failed")
os.makedirs(os.path.join(ROOT, "data", "history"), exist_ok=True)
out = os.path.join(ROOT, "data", "history", "spacetrack_gp_history.csv"); header = None; rows = 0
with open(out, "w", newline="") as fo:
    for k, b in enumerate(batches):
        if k: time.sleep(6)
        txt = op.open(url(b), timeout=180).read().decode("utf-8", "replace").splitlines()
        if not txt or "NORAD_CAT_ID" not in txt[0]: sys.exit(f"batch {k}: unexpected reply: {' '.join(txt[:2])[:200]}")
        if header is None: header = txt[0]; fo.write(header + "\n")
        elif txt[0] != header: sys.exit(f"batch {k}: column layout changed")
        for line in txt[1:]: fo.write(line + "\n"); rows += 1
        print(f"batch {k + 1}/{len(batches)}: {len(txt) - 1} rows")
try: op.open(BASE + "/ajaxauth/logout", timeout=30)
except Exception: pass
print(f"wrote {rows} element sets to {out}")
