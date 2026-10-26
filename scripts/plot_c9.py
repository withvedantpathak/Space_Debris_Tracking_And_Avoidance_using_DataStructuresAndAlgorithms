#!/usr/bin/env python3
"""C9 figures. Every plotted number is read from a file in results/ (or data/raw/ for V06) produced by a script in scripts/.
Schematics V01-V04 are drawings and carry no data. Usage: python3 scripts/plot_c9.py [V01 V13 ...]   (no argument = all)
Output: paper/figures/<id>_<name>.png (200 dpi) and .pdf.  Exports needed first: scripts/export_c9_data.sh"""
import sys, os, glob
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
os.chdir(ROOT)
R = "results/"
OUT = "paper/figures/"
os.makedirs(OUT, exist_ok=True)
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00", purple="#CC79A7", sky="#56B4E9", grey="#6b6b6b", dark="#222222")
METHOD_COL = dict(brute=C["grey"], kd=C["blue"], grid=C["orange"], alt=C["purple"])
METHOD_LBL = dict(brute="Brute force", kd="k-d tree", grid="Hashed grid", alt="Altitude filter")
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25,
                     "figure.dpi": 100, "savefig.bbox": "tight", "legend.frameon": False, "axes.titlesize": 11, "axes.titleweight": "bold"})
FIGS = {}
def fig(vid):
    def deco(f): FIGS[vid] = f; return f
    return deco
def save(f, vid, name):
    for ext, kw in (("png", dict(dpi=200)), ("pdf", {})):
        f.savefig(f"{OUT}{vid}_{name}.{ext}", **kw)
    plt.close(f); print("wrote", f"{OUT}{vid}_{name}.png")
def src(ax_or_fig, text):
    f = ax_or_fig.figure if hasattr(ax_or_fig, "figure") else ax_or_fig
    f.text(0.01, -0.05, "Source: " + text, fontsize=7, color=C["grey"], ha="left", va="top")

# ---------------------------------------------------------------- A. problem and design (schematics)
def box(ax, x, y, w, h, title, sub, col):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=col, ec=C["dark"], lw=1, alpha=0.95))
    ax.text(x + w / 2, y + h * 0.68, title, ha="center", va="center", fontsize=9.5, fontweight="bold", color="white")
    ax.text(x + w / 2, y + h * 0.30, sub, ha="center", va="center", fontsize=6.6, color="white", linespacing=1.4)

@fig("V01")
def v01():
    f, ax = plt.subplots(figsize=(11.6, 3.4)); ax.set_xlim(0, 11.6); ax.set_ylim(0, 3.4); ax.axis("off")
    steps = [("Catalogue", "OMM CSV by NORAD id\nhash table\ndynamic array", C["dark"]),
             ("Propagate", "SGP4 (real data)\nKepler + J2 (synthetic)\nstate at each step", C["grey"]),
             ("Index", "k-d tree or\nhashed grid\nradius query per step", C["blue"]),
             ("Refine", "bisection on g(t)=r.v\n(or linear / hybrid)\nTCA and miss distance", C["green"]),
             ("Graph", "CSR graph\nBFS/DFS: queue, stack\nunion-find", C["orange"]),
             ("Rank + mitigate", "heap Top-K\nAVL timeline\ngreedy vs exact cover", C["red"])]
    w, g = 1.62, 0.27
    for i, (t, s, c) in enumerate(steps):
        x = 0.15 + i * (w + g); box(ax, x, 1.1, w, 1.4, t, s, c)
        if i: ax.annotate("", xy=(x, 1.8), xytext=(x - g, 1.8), arrowprops=dict(arrowstyle="-|>", color=C["dark"], lw=1.4))
    ax.text(5.8, 3.05, "Conjunction-screening pipeline and the data structure behind each stage", ha="center", fontsize=11, fontweight="bold")
    ax.text(5.8, 0.55, "Exactness guarantee: a proven Lipschitz skip bound means the index never drops a pair that brute force would report (V03).",
            ha="center", fontsize=8.5, color=C["grey"])
    save(f, "V01", "pipeline")

@fig("V02")
def v02():
    rows = [("Dynamic array", "event lists, CSR graph", "push O(1) amortised"),
            ("Hash table", "NORAD-id lookup, grid cells, timeline", "get/put O(1) expected"),
            ("Binary heap", "Top-K closest, greedy max-degree", "push/pop O(log n)"),
            ("Queue", "BFS components, hop distance", "O(1)"),
            ("Stack", "iterative DFS, undo in exact search", "O(1)"),
            ("AVL tree", "time-ordered conjunction index", "insert/find O(log n)"),
            ("k-d tree", "3-D radius query for candidate pairs", "build O(n log n); query ~O(log n + k) typical"),
            ("Uniform grid (hash)", "cell lookup for candidate pairs", "build O(n); query O(1) cells + k typical"),
            ("Union-find (DSU)", "cluster cross-check", "near-constant per operation"),
            ("Greedy / branch and bound", "vertex-cover mitigation", "greedy O(m log n); exact exponential worst case")]
    f, ax = plt.subplots(figsize=(12, 4.3)); ax.axis("off")
    t = ax.table(cellText=rows, colLabels=["Structure / technique", "Used for", "Cost"], loc="center", cellLoc="left", colWidths=[0.22, 0.34, 0.44], bbox=[0, 0, 1, 0.92])
    t.auto_set_font_size(False); t.set_fontsize(8.5)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if r == 0: cell.set_facecolor(C["dark"]); cell.get_text().set_color("white"); cell.get_text().set_fontweight("bold")
        elif r % 2 == 0: cell.set_facecolor("#f3f3f3")
    ax.set_title("Hand-written data structures and where the pipeline uses them", loc="left")
    src(f, "docs/dsa_coverage.md (use), standard textbook costs (cost column)")
    save(f, "V02", "data_structure_map")

@fig("V03")
def v03():
    f, ax = plt.subplots(figsize=(8.2, 3.9)); h = 10.0; L = 0.25; thr = 5.0
    t = np.linspace(0, h, 100)
    ax.axhline(thr, color=C["red"], ls="--", lw=1.2); ax.text(0.15, thr + 0.12, "threshold", color=C["red"], va="bottom", fontsize=8)
    for d0, col, lab in ((9.5, C["blue"], "skipped"), (6.0, C["green"], "kept: bound reaches threshold")):
        lb = d0 - L * t
        ax.fill_between(t, lb, d0 + L * t, color=col, alpha=0.12)
        ax.plot(t, lb, color=col, lw=1.6, label=f"{lab}  (d0 = {d0:g})")
        ax.plot([0], [d0], "o", color=col)
    ax.annotate("d0 - L h: the lowest the\ndistance can reach", xy=(h, 9.5 - L * h), xytext=(6.4, 11.0), fontsize=8, color=C["blue"], arrowprops=dict(arrowstyle="->", color=C["blue"]))
    ax.set_xlim(0, h * 1.02); ax.set_ylim(2.5, 12.3); ax.set_xlabel("time within one step (0 to h)"); ax.set_ylabel("distance between the pair (arbitrary units)")
    ax.set_title("Why a pair can be skipped without risk: distance changes by at most L per unit time")
    ax.text(0.2, 3.0, "Skip the pair when  d0 > threshold + L h,  with L = vmax_i + vmax_j", fontsize=8.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=8); ax.text(0.0, -0.33, "Schematic, not data.", transform=ax.transAxes, fontsize=7, color=C["grey"])
    save(f, "V03", "lipschitz_skip")

@fig("V04")
def v04():
    rng = np.random.default_rng(7); P = rng.random((26, 2)) * 10
    f, (a, b) = plt.subplots(1, 2, figsize=(9.6, 4.3))
    def split(idx, x0, x1, y0, y1, d, depth):
        if len(idx) <= 2 or depth > 4: return
        ax_ = d % 2; v = np.sort(P[idx, ax_]); m = v[len(v) // 2]
        if ax_ == 0: a.plot([m, m], [y0, y1], color=C["blue"], lw=1)
        else: a.plot([x0, x1], [m, m], color=C["blue"], lw=1)
        lo = [i for i in idx if P[i, ax_] < m]; hi = [i for i in idx if P[i, ax_] >= m]
        if ax_ == 0: split(lo, x0, m, y0, y1, d + 1, depth + 1); split(hi, m, x1, y0, y1, d + 1, depth + 1)
        else: split(lo, x0, x1, y0, m, d + 1, depth + 1); split(hi, x0, x1, m, y1, d + 1, depth + 1)
    a.scatter(P[:, 0], P[:, 1], s=14, color=C["dark"], zorder=3); split(list(range(len(P))), 0, 10, 0, 10, 0, 0)
    q = (4.6, 5.2); a.add_patch(Circle(q, 1.9, fc="none", ec=C["red"], lw=1.6)); a.plot(*q, "x", color=C["red"])
    a.set_title("k-d tree: median splits", fontsize=10); a.set_xlim(0, 10); a.set_ylim(0, 10); a.set_aspect("equal")
    cs = 2.0
    for k in range(6): b.axhline(k * cs, color=C["orange"], lw=0.8); b.axvline(k * cs, color=C["orange"], lw=0.8)
    b.scatter(P[:, 0], P[:, 1], s=14, color=C["dark"], zorder=3); b.add_patch(Circle(q, 1.9, fc="none", ec=C["red"], lw=1.6)); b.plot(*q, "x", color=C["red"])
    for i in range(2, 4):
        for j in range(1, 4): b.add_patch(Rectangle((i * cs, j * cs), cs, cs, fc=C["orange"], alpha=0.12, lw=0))
    b.set_title("Hashed grid: look at 3x3 cells", fontsize=10); b.set_xlim(0, 10); b.set_ylim(0, 10); b.set_aspect("equal")
    f.text(0.5, -0.02, "Schematic with random 2-D points (the real index is 3-D).", ha="center", fontsize=7, color=C["grey"])
    save(f, "V04", "spatial_index_sketch")

# ---------------------------------------------------------------- B. data
def load_real():
    """Real catalogue: four CelesTrak files, one row per NORAD id (newest epoch). Returns DataFrame with source group, altitude, inclination."""
    names = [("cosmos-2251-debris", "Cosmos-2251 debris"), ("fengyun-1c-debris", "Fengyun-1C debris"), ("iridium-33-debris", "Iridium-33 debris"), ("active", "Active satellites and others")]
    fr = []
    for fn, lab in names:
        d = pd.read_csv(f"data/raw/{fn}.csv"); d["group"] = lab; fr.append(d)
    d = pd.concat(fr, ignore_index=True)
    d = d.sort_values(["NORAD_CAT_ID", "EPOCH"]).groupby("NORAD_CAT_ID", as_index=False).last()   # newest epoch wins (group = file of that row)
    mu, Re = 398600.4418, 6378.137
    n = d.MEAN_MOTION.values * 2 * np.pi / 86400
    a = (mu / n ** 2) ** (1 / 3)
    d["alt"] = a - Re; d["perigee"] = a * (1 - d.ECCENTRICITY) - Re; d["apogee"] = a * (1 + d.ECCENTRICITY) - Re
    return d

@fig("V05")
def v05():
    d = pd.read_csv(R + "dataset_stats.csv")
    leo, geo, hi = d.leo_share.values, d.geo_share.values, d.high_ecc_share.values
    oth = 1 - leo - geo - hi
    labels = ["Real (CelesTrak)", "Synthetic (20,000)"]
    f, ax = plt.subplots(figsize=(7.5, 2.9)); left = np.zeros(2)
    for vals, lab, col in ((leo, "LEO", C["blue"]), (geo, "GEO", C["orange"]), (hi, "High eccentricity", C["red"]), (oth, "MEO / other", C["grey"])):
        ax.barh(labels, vals * 100, left=left * 100, color=col, label=lab)
        for i, v in enumerate(vals):
            if v > 0.02: ax.text((left[i] + v / 2) * 100, i, f"{v*100:.1f}%", ha="center", va="center", color="white", fontsize=8)
        left += vals
    ax.set_xlabel("share of objects (%)"); ax.set_xlim(0, 100); ax.invert_yaxis(); ax.grid(False)
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.28), fontsize=8)
    ax.set_title("Orbit-class mix: real and synthetic catalogues")
    src(f, "results/dataset_stats.csv"); save(f, "V05", "catalogue_composition")

@fig("V06")
def v06():
    d = load_real(); m = (d.perigee >= 175) & (d.apogee <= 2000)
    f, ax = plt.subplots(figsize=(8.4, 4.8))
    order = ["Active satellites and others", "Cosmos-2251 debris", "Fengyun-1C debris", "Iridium-33 debris"]
    cols = [C["grey"], C["orange"], C["red"], C["blue"]]
    for g, c in zip(order, cols):
        s = d[m & (d.group == g)]; ax.scatter(s.INCLINATION, s.alt, s=3 if g.startswith("Active") else 5, color=c, alpha=0.5 if g.startswith("Active") else 0.8, label=f"{g} ({len(s):,})", lw=0)
    ax.set_xlabel("inclination (deg)"); ax.set_ylabel("mean altitude (km)"); ax.set_ylim(150, 2050)
    ax.set_title(f"LEO-scope catalogue: {int(m.sum()):,} objects, altitude vs inclination")
    ax.legend(markerscale=4, fontsize=8, loc="upper right")
    src(f, "data/raw/*.csv (mean altitude from mean motion; newest element set per NORAD id; scope perigee >= 175 km, apogee <= 2000 km)"); save(f, "V06", "altitude_inclination")

@fig("V07")
def v07():
    p = pd.read_csv(R + "c9_positions.csv"); meta = pd.read_csv(R + "c9_pair_meta.csv", header=None)
    A, B = int(meta.iloc[0, 1]), int(meta.iloc[0, 2])
    r = np.sqrt(p.x_km ** 2 + p.y_km ** 2 + p.z_km ** 2); p = p[(r > 6500) & (r < 8400)].reset_index(drop=True)
    az, el = np.radians(35), np.radians(28)                      # fixed orthographic view; depth > 0 is towards the viewer
    def proj(x, y, z):
        x1 = x * np.cos(az) - y * np.sin(az); y1 = x * np.sin(az) + y * np.cos(az)
        return x1, z * np.cos(el) - y1 * np.sin(el), y1 * np.cos(el) + z * np.sin(el)
    X, Y, D = proj(p.x_km.values, p.y_km.values, p.z_km.values); Re = 6378.137
    hidden = (D < 0) & (np.hypot(X, Y) < Re)                      # behind the Earth disk
    f, ax = plt.subplots(figsize=(7.4, 7.2)); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(Circle((0, 0), Re, fc="#bcd0e4", ec="#7f9bb8", lw=1, zorder=1))
    ax.scatter(X[~hidden & (D < 0)], Y[~hidden & (D < 0)], s=0.7, color=C["grey"], alpha=0.5, lw=0, zorder=0)
    ax.scatter(X[D >= 0], Y[D >= 0], s=0.9, color=C["dark"], alpha=0.55, lw=0, zorder=2)
    sel = p.norad_id.isin([A, B]).values
    ax.scatter(X[sel], Y[sel], s=160, facecolor="none", edgecolor=C["red"], lw=2, zorder=3)
    ax.annotate(f"closest crossing pair\n{A} and {B}\n(miss {float(meta.iloc[0, 6]):.2f} km)", xy=(X[sel].mean(), Y[sel].mean()), xytext=(Re * 0.95, -Re * 1.25), color=C["red"], fontsize=9, arrowprops=dict(arrowstyle="->", color=C["red"]), zorder=4)
    ax.set_xlim(-8800, 8800); ax.set_ylim(-8800, 8800)
    ax.set_title(f"{len(p):,} LEO objects at the moment of the closest crossing encounter (1 h run)", fontsize=10)
    src(f, "results/c9_positions.csv (dlscms prop, SGP4), results/c9_pair_meta.csv; fixed orthographic view, Earth drawn to scale"); save(f, "V07", "orbit_view")

# ---------------------------------------------------------------- C. correctness
@fig("V08")
def v08():
    d = pd.read_csv(R + "propagation_error.csv"); d = d[d.orbit_class == "LEO"]
    f, ax = plt.subplots(figsize=(7.4, 4.3))
    for m, c in (("two-body", C["grey"]), ("two-body + J2", C["orange"]), ("two-body + J2 + drag", C["blue"])):
        s = d[d.model == m].sort_values("dt_h"); ax.plot(s.dt_h, s.median_km, "o-", color=c, label=f"{m} (median)")
        ax.plot(s.dt_h, s.p95_km, "--", color=c, alpha=0.6, lw=1)
    ax.set_yscale("log"); ax.set_xlabel("time from epoch (h)"); ax.set_ylabel("position error vs SGP4 (km)")
    ax.set_title("Simplified model against SGP4 for LEO objects (dashed = 95th percentile)"); ax.legend(fontsize=8)
    src(f, "results/propagation_error.csv"); save(f, "V08", "propagation_error")

@fig("V09")
def v09():
    d = pd.read_csv(R + "sgp4_validation.csv")
    f, ax = plt.subplots(figsize=(6.6, 3.8)); x = np.arange(len(d))
    ax.bar(x, d.max_pos_diff_km * 1e6, color=C["blue"]); ax.set_xticks(x); ax.set_xticklabels(d.t)
    ax.axhline(1000, color=C["red"], ls="--"); ax.text(len(d) - 0.5, 1150, "1 mm", color=C["red"], ha="right", fontsize=8)
    ax.set_yscale("log"); ax.set_ylim(0.1, 3000); ax.set_ylabel("max position difference (mm)"); ax.set_xlabel("time after epoch")
    ax.set_title(f"Hand-written SGP4 vs reference implementation ({int(d.compared.iloc[0]):,} objects)")
    src(f, "results/sgp4_validation.csv"); save(f, "V09", "sgp4_validation")

@fig("V10")
def v10():
    d = pd.read_csv(R + "event_agreement.csv"); d = d[~d.event_class.str.contains("precision")]
    f, axs = plt.subplots(1, 2, figsize=(9.6, 3.9), sharey=True)
    for ax, thr in zip(axs, (5, 10)):
        for cl, c in (("all", C["dark"]), ("crossing (v_rel > 1 km/s)", C["blue"]), ("slow / co-orbiting (v_rel <= 1 km/s)", C["orange"])):
            s = d[(d.event_class == cl) & (d.sgp4_threshold_km == thr)].sort_values("our_threshold_km")
            ax.plot(s.our_threshold_km, s.recall * 100, "o-", color=c, label=cl)
        ax.set_title(f"SGP4 events below {thr} km"); ax.set_xlabel("screening distance used (km)")
    axs[0].set_ylabel("share of SGP4 events found (%)"); axs[0].legend(fontsize=7.5, loc="lower right")
    src(f, "results/event_agreement.csv"); save(f, "V10", "event_agreement")

@fig("V11")
def v11():
    c = pd.read_csv(R + "c9_pair_curve.csv"); m = pd.read_csv(R + "c9_pair_meta.csv", header=None); tca, miss = float(m.iloc[0, 4]), float(m.iloc[0, 6])
    r = c[["dx", "dy", "dz"]].values; v = c[["dvx", "dvy", "dvz"]].values
    dist = np.linalg.norm(r, axis=1); g = (r * v).sum(1)
    f, (a, b) = plt.subplots(2, 1, figsize=(7.4, 5.2), sharex=True)
    a.plot(c.dt_s, dist, color=C["blue"]); a.set_ylabel("distance (km)"); a.set_title("Closest-approach refinement for the closest crossing pair of the 1 h run")
    b.plot(c.dt_s, g, color=C["green"]); b.axhline(0, color=C["dark"], lw=0.8); b.set_ylabel("g(t) = r . v  (km^2/s)"); b.set_xlabel("time since window start (s)")
    for ax in (a, b):
        ax.axvline(tca, color=C["red"], ls="--", lw=1)
        for t in (1620, 1630):
            if c.dt_s.min() <= t <= c.dt_s.max(): ax.axvline(t, color=C["grey"], lw=0.8, ls=":")
    a.annotate(f"TCA {tca:.3f} s\nmiss {miss:.4f} km", xy=(tca, miss), xytext=(tca + 6, dist.max() * 0.55), arrowprops=dict(arrowstyle="->"), fontsize=8, color=C["red"])
    b.text(1631, -g.max() * 0.75, "10 s screening steps (dotted):\nsign change of g(t) inside one step\n=> bisect for the root", fontsize=8)
    src(f, "results/c9_pair_curve.csv, c9_pair_meta.csv (dlscms prop), c6_1h_topk_crossing.csv"); save(f, "V11", "refinement_curve")

@fig("V12")
def v12():
    s = pd.read_csv(R + "c7_scaling.csv"); b5 = pd.read_csv(R + "bench_c5.csv")
    rows = []
    for n in sorted(s.n.unique()):
        row = [f"{n:,} (synthetic)"]
        for m in ("kd", "grid"):
            q = s[(s.n == n) & (s.method == m)]; row.append(f"{(q.identical_to_bruteforce == 'yes').sum()} / {len(q)}")
        rows.append(row)
    for n, lab in ((20000, "20,000 (synthetic)"), (18191, "18,191 (real, LEO scope)")):
        row = [lab]
        for m in ("kd", "grid"):
            q = b5[(b5.n == n) & (b5.method == m)]; row.append(f"{(q.identical_to_bruteforce == 'yes').sum()} / {len(q)}")
        rows.append(row)
    f, ax = plt.subplots(figsize=(6.2, 3.2)); ax.axis("off")
    t = ax.table(cellText=rows, colLabels=["Catalogue size N", "k-d tree: runs identical", "Grid: runs identical"], loc="center", cellLoc="center")
    t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.5)
    for (r_, c_), cell in t.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if r_ == 0: cell.set_facecolor(C["dark"]); cell.get_text().set_color("white")
    ax.set_title("Fast methods produce the same event file as brute force (byte-identical)", fontsize=10)
    src(f, "results/c7_scaling.csv, results/bench_c5.csv. Runs = repeats x seeds (brute force once for N >= 10,000; one run at 18,191 and 20,000)"); save(f, "V12", "equality_with_bruteforce")

# ---------------------------------------------------------------- D. performance
def scaling_table():
    """Mean/min/max CPU, pair checks and RSS per (n, method): C7 scaling study + 50k point + 20k synthetic and real from bench_c5."""
    s = pd.read_csv(R + "c7_scaling.csv").assign(src="c7")
    s5 = pd.read_csv(R + "c7_scaling_50k.csv").assign(src="c7_50k")
    b = pd.read_csv(R + "bench_c5.csv"); b = b[(b.n == 20000) | (b.n == 18191)]
    b = b.rename(columns={"cpu_seconds": "cpu_seconds"}).assign(src="bench_c5", peak_rss_kb=np.nan)[["n", "method", "cpu_seconds", "pair_checks", "peak_rss_kb", "src"]]
    a = pd.concat([s[["n", "method", "cpu_seconds", "pair_checks", "peak_rss_kb", "src"]], s5[["n", "method", "cpu_seconds", "pair_checks", "peak_rss_kb", "src"]], b], ignore_index=True)
    g = a.groupby(["n", "method"]).agg(cpu=("cpu_seconds", "mean"), lo=("cpu_seconds", "min"), hi=("cpu_seconds", "max"), pc=("pair_checks", "mean"), rss=("peak_rss_kb", "mean")).reset_index()
    return g

@fig("V13")
def v13():
    g = scaling_table(); f, ax = plt.subplots(figsize=(7.4, 4.8))
    for m in ("brute", "kd", "grid"):
        s = g[(g.method == m) & (g.n != 18191)].sort_values("n")
        ax.errorbar(s.n, s.cpu, yerr=[s.cpu - s.lo, s.hi - s.cpu], fmt="o-", color=METHOD_COL[m], label=METHOD_LBL[m], capsize=2)
    real = g[g.n == 18191]
    for _, r_ in real.iterrows(): ax.plot(r_.n, r_.cpu, "s", mfc="none", color=METHOD_COL[r_.method], ms=8)
    ax.plot([], [], "s", mfc="none", color=C["dark"], label="real catalogue (LEO scope, 18,191)")
    s1 = g[(g.method == "brute") & (g.n == 1000)].cpu.iloc[0]; xs = np.array([1000, 20000]); ax.plot(xs, s1 * (xs / 1000) ** 2, ":", color=C["grey"], lw=1); ax.text(7000, s1 * 28, "N^2 guide", color=C["grey"], fontsize=8, rotation=32)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("catalogue size N (objects)"); ax.set_ylabel("CPU time, 1 h window (s)")
    ax.set_title("Run time against catalogue size (synthetic; bars = min to max of repeats and seeds)"); ax.legend(fontsize=8)
    ax.text(0.01, 0.01, "No brute-force run at N = 50,000 (see L51).", transform=ax.transAxes, fontsize=7.5, color=C["grey"])
    src(f, "results/c7_scaling.csv, c7_scaling_50k.csv, bench_c5.csv"); save(f, "V13", "runtime_vs_n")

@fig("V14")
def v14():
    g = scaling_table(); f, ax = plt.subplots(figsize=(7.4, 4.6))
    for m in ("brute", "kd", "grid"):
        s = g[(g.method == m) & (g.n != 18191)].sort_values("n"); ax.plot(s.n, s.pc, "o-", color=METHOD_COL[m], label=METHOD_LBL[m])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("catalogue size N (objects)"); ax.set_ylabel("pair distance checks, 1 h window")
    n20 = g[g.n == 20000].set_index("method").pc; ratio = n20["brute"] / n20["kd"]
    ax.annotate(f"{ratio:,.0f}x fewer checks at N = 20,000", xy=(20000, n20["kd"]), xytext=(1400, n20["kd"] * 0.05), fontsize=9, arrowprops=dict(arrowstyle="->"))
    ax.set_title("Pairs examined: same events with far fewer checks"); ax.legend(fontsize=8)
    src(f, "results/c7_scaling.csv, c7_scaling_50k.csv, bench_c5.csv"); save(f, "V14", "pair_checks_vs_n")

@fig("V15")
def v15():
    g = scaling_table(); g = g.dropna(subset=["rss"]); c = pd.read_csv(R + "c7_cell.csv")
    f, (a, b) = plt.subplots(1, 2, figsize=(9.6, 4), gridspec_kw=dict(width_ratios=[1.3, 1]))
    for m in ("brute", "kd", "grid"):
        s = g[g.method == m].sort_values("n"); a.plot(s.n, s.rss / 1024, "o-", color=METHOD_COL[m], label=METHOD_LBL[m])
    a.set_xscale("log"); a.set_yscale("log"); a.set_xticks([1000, 2000, 5000, 10000, 50000]); a.set_xticklabels(["1k", "2k", "5k", "10k", "50k"]); a.set_yticks([4, 6, 10, 20, 40, 60]); a.set_yticklabels(["4", "6", "10", "20", "40", "60"]); a.minorticks_off(); a.set_xlabel("catalogue size N"); a.set_ylabel("peak resident memory (MB)"); a.legend(fontsize=8); a.set_title("Peak memory against N")
    gr = c[c.method == "grid"]; b.plot(gr.cell_scale, gr.peak_rss_kb / 1024, "o-", color=C["orange"], label="grid"); b.axhline(c[c.method == "kd"].peak_rss_kb.iloc[0] / 1024, color=C["blue"], ls="--", label="k-d tree")
    b.set_xscale("log"); b.set_ylim(0, 30); b.set_xlabel("grid cell scale"); b.set_ylabel("peak memory (MB)"); b.legend(fontsize=8); b.set_title("Memory against cell size (real data, 10 min)")
    src(f, "results/c7_scaling.csv, c7_scaling_50k.csv, c7_cell.csv"); save(f, "V15", "memory")

@fig("V16")
def v16():
    c = pd.read_csv(R + "c7_cell.csv"); gr = c[c.method == "grid"]; kd = c[c.method == "kd"].cpu_seconds.iloc[0]
    f, ax = plt.subplots(figsize=(6.6, 4)); ax.plot(gr.cell_scale, gr.cpu_seconds, "o-", color=C["orange"], label="grid")
    ax.axhline(kd, color=C["blue"], ls="--", label="k-d tree"); b = gr.loc[gr.cpu_seconds.idxmin()]; ax.plot(b.cell_scale, b.cpu_seconds, "o", ms=11, mfc="none", color=C["red"])
    ax.annotate(f"best: scale {b.cell_scale:g}\n{b.cpu_seconds:.2f} s", xy=(b.cell_scale, b.cpu_seconds), xytext=(b.cell_scale * 1.3, b.cpu_seconds * 2.2), arrowprops=dict(arrowstyle="->"), fontsize=8)
    ax.set_xscale("log"); ax.set_xlabel("cell size relative to the query reach"); ax.set_ylabel("CPU time (s)"); ax.legend(fontsize=8)
    ax.set_title("Grid cell size: identical events and pair checks at every scale, time differs"); src(f, "results/c7_cell.csv"); save(f, "V16", "cell_size")

@fig("V17")
def v17():
    s = pd.read_csv(R + "c7_step.csv"); f, (a, b) = plt.subplots(1, 2, figsize=(9.2, 3.9))
    a.plot(s.step_s, s.events, "o-", color=C["green"]); a.set_ylim(0, s.events.max() * 1.3); a.set_xlabel("screening step (s)"); a.set_ylabel("events found"); a.set_title("Events: unchanged from 5 s to 60 s")
    b.plot(s.step_s, s.cpu_seconds, "o-", color=C["red"], label="CPU time (s)"); b.set_xlabel("screening step (s)"); b.set_ylabel("CPU time (s)")
    b2 = b.twinx(); b2.plot(s.step_s, s.candidates / 1e6, "s--", color=C["grey"], label="candidate intervals (millions)"); b2.set_ylabel("candidate intervals (millions)"); b2.grid(False)
    b.set_title("A larger step only adds cost"); h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels(); b.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8)
    src(f, "results/c7_step.csv"); save(f, "V17", "step_study")

@fig("V18")
def v18():
    b = pd.read_csv(R + "bench_c5.csv"); b = b[b.n != 18191]
    f, (a, c) = plt.subplots(1, 2, figsize=(9.8, 4))
    for m in ("brute", "kd", "grid", "alt"):
        s = b[b.method == m].sort_values("n"); a.plot(s.n, s.cpu_seconds, "o-", color=METHOD_COL[m], label=METHOD_LBL[m]); c.plot(s.n, s.pair_checks, "o-", color=METHOD_COL[m], label=METHOD_LBL[m])
    for ax, yl in ((a, "CPU time (s)"), (c, "pair distance checks")):
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("catalogue size N"); ax.set_ylabel(yl)
    a.legend(fontsize=8); a.set_title("Altitude filter is slower than brute force"); c.set_title("...and removes only a few pairs")
    src(f, "results/bench_c5.csv"); save(f, "V18", "altitude_filter")

@fig("V19")
def v19():
    b = pd.read_csv(R + "bench_c1.csv"); g = b.groupby(["method", "n"]).agg(build=("build_ms", "mean"), query=("query_ms", "mean")).reset_index()
    col = dict(brute=C["grey"], grid=C["orange"], kdtree=C["blue"]); lab = dict(brute="Brute force", grid="Hashed grid", kdtree="k-d tree")
    f, (a, c) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    for m in ("brute", "grid", "kdtree"):
        s = g[g.method == m].sort_values("n"); a.plot(s.n, s["query"], "o-", color=col[m], label=lab[m])
        if m != "brute": c.plot(s.n, s.build, "o-", color=col[m], label=lab[m])
    for ax, yl in ((a, "time to find all pairs within 50 km (ms)"), (c, "index build time (ms)")): ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("number of points N"); ax.set_ylabel(yl)
    a.legend(fontsize=8); a.set_title("Fixed-radius pair search (mean of 3 seeds)"); c.set_title("Index construction")
    src(f, "results/bench_c1.csv"); save(f, "V19", "structure_microbenchmarks")

# ---------------------------------------------------------------- E. graph, clusters, mitigation
def components_of(edges):
    par = {}
    def find(x):
        while par.setdefault(x, x) != x: par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in edges: par[find(a)] = find(b)
    comp = {}
    for x in list(par): comp.setdefault(find(x), []).append(x)
    return sorted(comp.values(), key=len, reverse=True)

@fig("V20")
def v20():
    ev = pd.read_csv(R + "c7_events_thr5.csv"); ev["a"] = ev[["norad_id_a", "norad_id_b"]].min(axis=1); ev["b"] = ev[["norad_id_a", "norad_id_b"]].max(axis=1)
    pr = ev.groupby(["a", "b"]).agg(rel=("rel_speed_kms", "max"), miss=("miss_km", "min")).reset_index()
    big = components_of(zip(pr.a, pr.b))[0]; nodes = sorted(big); idx = {n: i for i, n in enumerate(nodes)}
    E = pr[pr.a.isin(big) & pr.b.isin(big)]; n = len(nodes); deg = np.zeros(n)
    for _, r_ in E.iterrows(): deg[idx[r_.a]] += 1; deg[idx[r_.b]] += 1
    rng = np.random.default_rng(3); pos = rng.random((n, 2)); k = np.sqrt(1.0 / n); t = 0.1
    for it in range(400):                       # Fruchterman-Reingold layout
        d = pos[:, None, :] - pos[None, :, :]; dist = np.linalg.norm(d, axis=2) + 1e-9
        disp = (d / dist[..., None] * (k * k / dist)[..., None]).sum(1)
        for _, r_ in E.iterrows():
            i, j = idx[r_.a], idx[r_.b]; dd = pos[i] - pos[j]; L_ = np.linalg.norm(dd) + 1e-9; fvec = dd / L_ * (L_ * L_ / k)
            disp[i] -= fvec; disp[j] += fvec
        ln = np.linalg.norm(disp, axis=1)[:, None] + 1e-9; pos += disp / ln * np.minimum(ln, t); t *= 0.985
    cov = pd.read_csv(R + "c7_thr5_cover.csv").set_index("norad_id"); best = [bool(cov.in_best_cover.get(x, 0)) for x in nodes]
    f, ax = plt.subplots(figsize=(7.2, 6.4))
    for _, r_ in E.iterrows():
        p1, p2 = pos[idx[r_.a]], pos[idx[r_.b]]; crossing = r_.rel >= 1.0
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=C["red"] if crossing else C["grey"], lw=1.4 if crossing else 0.9, alpha=0.8, zorder=1)
    sc = ax.scatter(pos[:, 0], pos[:, 1], s=60 + 30 * deg, c=deg, cmap="viridis", edgecolor=[C["red"] if b_ else "white" for b_ in best], linewidths=[2.2 if b_ else 0.8 for b_ in best], zorder=2)
    f.colorbar(sc, ax=ax, label="number of conjunction partners (degree)", shrink=0.7)
    ax.plot([], [], color=C["red"], label="crossing pair (v_rel >= 1 km/s)"); ax.plot([], [], color=C["grey"], label="slow / co-orbiting pair"); ax.plot([], [], "o", mfc="none", mec=C["red"], mew=2, label="in best (exact) mitigation set")
    ax.legend(fontsize=8, loc="lower left"); ax.set_axis_off(); ax.grid(False)
    ax.set_title(f"Largest conjunction cluster at 5 km: {n} objects, {len(E)} pairs (1 h window)")
    src(f, "results/c7_events_thr5.csv, c7_thr5_cover.csv"); save(f, "V20", "cluster_graph")

@fig("V21")
def v21():
    c = pd.read_csv(R + "c7_thr10_components.csv"); t = pd.read_csv(R + "c7_threshold.csv")
    f, (a, b) = plt.subplots(1, 2, figsize=(9.8, 4))
    bins = np.unique(np.round(np.logspace(0, np.log10(c.nodes.max() + 1), 18)).astype(int)); h, e = np.histogram(c.nodes, bins=bins)
    a.bar(e[:-1], h, width=np.diff(e), align="edge", color=C["blue"], edgecolor="white"); a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("objects in cluster"); a.set_ylabel("number of clusters")
    a.set_title(f"Cluster sizes at 10 km (largest {int(c.nodes.max()):,})", fontsize=10)
    b.plot(t.threshold_km, t.largest_cluster, "o-", color=C["red"], label="largest cluster (objects)"); b.plot(t.threshold_km, t.components, "s--", color=C["blue"], label="number of clusters"); b.plot(t.threshold_km, t.objects, "^:", color=C["grey"], label="objects in any conjunction")
    b.set_yscale("log"); b.set_xlabel("screening distance (km)"); b.set_xticks(t.threshold_km); b.legend(fontsize=8); b.set_title("Giant cluster appears as distance grows", fontsize=10)
    src(f, "results/c7_thr10_components.csv, c7_threshold.csv"); save(f, "V21", "clusters_vs_threshold")

@fig("V22")
def v22():
    f, axs = plt.subplots(2, 1, figsize=(8.2, 5.6), sharey=False); f.subplots_adjust(hspace=0.5)
    for ax, h in zip(axs, ("1h", "6h")):
        t = pd.read_csv(R + f"c6_{h}_timeline.csv"); ax.bar(t.minutes_from_start / 60 + 5 / 60, t.conjunction_pairs_with_tca_in_bin, width=10 / 60, color=C["blue"])
        ax.set_ylabel("pairs per 10 min"); ax.set_title(f"{h} window, 10 km", loc="left", fontsize=9)
    axs[1].set_xlabel("hours from window start"); axs[0].set_title("Conjunction timeline (AVL time index): events per 10-minute bin\n1 h window, 10 km", loc="left", fontsize=10)
    axs[1].set_title("6 h window, 10 km", loc="left", fontsize=9)
    src(f, "results/c6_1h_timeline.csv, c6_6h_timeline.csv"); save(f, "V22", "timeline")

@fig("V23")
def v23():
    t = pd.read_csv(R + "c6_1h_topk_crossing.csv").head(10)
    nm = {}
    for fn in ("active", "cosmos-2251-debris", "fengyun-1c-debris", "iridium-33-debris"):
        d = pd.read_csv(f"data/raw/{fn}.csv", usecols=["OBJECT_NAME", "NORAD_CAT_ID"]); nm.update(dict(zip(d.NORAD_CAT_ID, d.OBJECT_NAME)))
    rows = [[int(r_["rank"]), f"{nm.get(int(r_.norad_id_a), '?')[:22]} ({int(r_.norad_id_a)})", f"{nm.get(int(r_.norad_id_b), '?')[:22]} ({int(r_.norad_id_b)})", f"{r_.tca_s/60:.1f}", f"{r_.miss_km:.3f}", f"{r_.rel_speed_kms:.1f}", int(r_.degree_a), int(r_.degree_b)] for _, r_ in t.iterrows()]
    f, ax = plt.subplots(figsize=(11, 3.9)); ax.axis("off")
    tb = ax.table(cellText=rows, colLabels=["#", "Object A", "Object B", "TCA (min)", "Miss (km)", "v_rel (km/s)", "deg A", "deg B"], loc="center", cellLoc="left", colWidths=[0.04, 0.27, 0.27, 0.09, 0.09, 0.1, 0.07, 0.07], bbox=[0, 0, 1, 0.92])
    tb.auto_set_font_size(False); tb.set_fontsize(8)
    for (r_, c_), cell in tb.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if r_ == 0: cell.set_facecolor(C["dark"]); cell.get_text().set_color("white")
    ax.set_title("Ten closest crossing encounters in the 1 h run (10 km screening, v_rel > 1 km/s)", loc="left")
    src(f, "results/c6_1h_topk_crossing.csv; names from data/raw/*.csv"); save(f, "V23", "top_k_crossing")

@fig("V24")
def v24():
    s = pd.read_csv(R + "c6_summary.csv"); s = s[s.variant == "reductions"].reset_index(drop=True); wins = ["1 h", "3 h", "6 h"]
    milp = [pd.read_csv(R + f"c6_{h}_verify.csv").milp_optimum.iloc[0] for h in ("1h", "3h", "6h")]
    f, ax = plt.subplots(figsize=(8, 4.2)); x = np.arange(3); w = 0.2
    series = [("Lower bound", s.lower_bound, C["grey"]), ("Independent optimum (MILP)", milp, C["green"]), ("Exact search", s.exact_cover, C["blue"]), ("Greedy", s.greedy_cover, C["orange"])]
    for i, (lab, v, c) in enumerate(series):
        ax.bar(x + (i - 1.5) * w, v, w, label=lab, color=c)
        if lab == "Greedy":
            for xi, (gv, mv) in enumerate(zip(v, milp)): ax.text(xi + 1.5 * w, gv + 50, f"+{(gv/mv-1)*100:.1f}%", ha="center", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels([f"{w_} window" for w_ in wins]); ax.set_ylabel("objects to manoeuvre"); ax.set_ylim(0, max(s.greedy_cover) * 1.15)
    ax.set_title("Mitigation: how many objects must move to clear all conjunctions (10 km)"); ax.legend(fontsize=8, loc="upper left")
    un = [(w_, int(a_ - b_), int(a_)) for w_, a_, b_ in zip(wins, s.nontrivial_components, s.exact_finished) if a_ > b_]
    f.text(0.5, -0.02, "Exact search stopped at its node limit in " + "; ".join(f"{k} of {t} clusters ({w_})" for w_, k, t in un) + ". Cover sizes in objects.", ha="center", fontsize=7.5, color=C["grey"])
    src(f, "results/c6_summary.csv, c6_*h_verify.csv"); save(f, "V24", "mitigation")

@fig("V25")
def v25():
    f, (a, b) = plt.subplots(1, 2, figsize=(10, 4)); fr = []
    for h in ("1h", "3h", "6h"): d = pd.read_csv(R + f"c6_{h}_components.csv"); d["win"] = h; fr.append(d)
    d = pd.concat(fr); d = d[d.nodes > 3]
    ok, no = d[d.exact_complete == 1], d[d.exact_complete == 0]
    a.scatter(ok.nodes, ok.exact_us / 1e3, s=10, color=C["blue"], alpha=0.6, label="exact search finished"); a.scatter(no.nodes, no.exact_us / 1e3, s=30, color=C["red"], marker="x", label="stopped at node limit")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("objects in cluster"); a.set_ylabel("exact search time (ms)"); a.legend(fontsize=8); a.set_title("Exact search effort per cluster (1, 3, 6 h)")
    s = pd.read_csv(R + "c6_summary.csv"); wins = ["1 h", "3 h", "6 h"]; x = np.arange(3)
    for i, (v, lab, c) in enumerate((("reductions", "with reductions", C["blue"]), ("branch_and_bound_only", "branch and bound only", C["orange"]))):
        q = s[s.variant == v].reset_index(drop=True); b.bar(x[:len(q)] + (i - 0.5) * 0.35, q.exact_seconds, 0.35, label=lab, color=c)
    b.set_xticks(x); b.set_xticklabels(wins); b.set_ylabel("exact search time (s)"); b.legend(fontsize=8); b.set_title("Effect of the reduction rules (6 h: branch and bound only not run)", fontsize=9)
    src(f, "results/c6_*h_components.csv, c6_summary.csv"); save(f, "V25", "exact_search_effort")

# ---------------------------------------------------------------- F. open questions and quality
@fig("V26")
def v26():
    l = pd.read_csv(R + "c8_l33_interp.csv"); st = pd.read_csv(R + "c7_step.csv"); ht = pd.read_csv(R + "c8_hybrid_timing.csv")
    f, (a, b, c) = plt.subplots(1, 3, figsize=(13, 3.9))
    for (sc, ls), meth in ((("all", "-"), None),):
        pass
    for meth, col in (("linear", C["red"]), ("hybrid", C["green"])):
        for sc, ls in (("all", "-"), ("crossing", "--")):
            q = l[(l.method == meth) & (l.scope == sc)].sort_values("step_s"); a.plot(q.step_s, q.f1 * 100, "o" + ls, color=col, label=f"{meth}, {sc} events")
    a.set_ylim(94.5, 100.3); a.set_xlabel("screening step (s)"); a.set_ylabel("F1 vs bisection reference (%)"); a.legend(fontsize=7.5, loc="center right"); a.set_title("Agreement with bisection")
    for meth, col in (("linear", C["red"]), ("hybrid", C["green"])):
        q = l[(l.method == meth) & (l.scope == "all")].sort_values("step_s"); b.plot(q.step_s, q.extra, "o-", color=col, label=meth)
    b.set_xlabel("screening step (s)"); b.set_ylabel("extra events (not in reference)"); b.legend(fontsize=8); b.set_title("Pure linear invents events for slow pairs")
    c.plot(st.step_s, st.cpu_seconds, "o-", color=C["grey"], label="bisection (default)")
    for meth, col in (("linear", C["red"]), ("hybrid", C["green"])):
        q = ht[ht.method == meth].sort_values("step_s"); c.plot(q.step_s, q.cpu_seconds, "o-", color=col, label=meth)
    c.set_yscale("log"); c.set_xlabel("screening step (s)"); c.set_ylabel("CPU time (s)"); c.legend(fontsize=8); c.set_title("Run time (real catalogue, 1 h, 5 km)")
    src(f, "results/c8_l33_interp.csv, c7_step.csv, c8_hybrid_timing.csv"); save(f, "V26", "linear_interpolation")

@fig("V27")
def v27():
    d = pd.read_csv(R + "c8_age_calibration.csv"); a_ = pd.read_csv(R + "c8_age_distribution.csv"); f, (a, b) = plt.subplots(1, 2, figsize=(11, 3.9)); f.subplots_adjust(wspace=0.5)
    a.plot(d.growth_km_per_day, d.largest_cluster, "o-", color=C["red"], label="largest cluster (objects)"); a.set_yscale("log"); a.set_xlabel("age growth g (km per day of element age, per object)"); a.set_ylabel("largest cluster")
    a2 = a.twinx(); a2.plot(d.growth_km_per_day, d.events, "s--", color=C["blue"], label="events"); a2.set_ylabel("events"); a2.grid(False)
    a.axvline(1.0, color=C["grey"], ls=":"); a.text(1.05, d.largest_cluster.min() * 1.05, "about 1 km/day: implied by the\nabstract of [5], not calibrated here", fontsize=7.5, va="bottom", color=C["grey"])
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); a.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left"); a.set_title("Sensitivity to the age-aware distance")
    p = a_[a_.statistic.str.startswith("age_percentile")]; lab = [s.split("_")[2] for s in p.statistic]
    b.bar(["P" + l for l in lab], p.value, color=C["green"]); b.bar(["max"], a_[a_.statistic == "age_max_days"].value, color=C["grey"])
    b.set_ylabel("element-set age (days)"); b.set_title("Most objects are fresh: age distribution"); 
    src(f, "results/c8_age_calibration.csv, c8_age_distribution.csv"); save(f, "V27", "age_aware_distance")

@fig("V28")
def v28():
    c = pd.read_csv(R + "c8_coverage.csv"); m = pd.read_csv(R + "c8_mutation.csv"); cc = c[c.file != "TOTAL"].sort_values("branch_pct")
    f, axs = plt.subplots(1, 3 if os.path.exists(R + "c8_test_counts.csv") else 2, figsize=(14 if os.path.exists(R + "c8_test_counts.csv") else 10.5, 4.8), gridspec_kw=dict(width_ratios=([1.7, 0.7, 1] if os.path.exists(R + "c8_test_counts.csv") else [2, 0.8])))
    a, b = axs[0], axs[1]; y = np.arange(len(cc)); a.barh(y - 0.2, cc.line_pct, 0.4, color=C["blue"], label="lines"); a.barh(y + 0.2, cc.branch_pct, 0.4, color=C["orange"], label="branches")
    a.set_yticks(y); a.set_yticklabels([s.replace("src/", "") for s in cc.file], fontsize=7.5); a.set_xlim(60, 101); a.set_xlabel("coverage (%)"); a.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)
    tot = c[c.file == "TOTAL"].iloc[0]; a.set_title(f"gcov coverage: {tot.line_pct:.1f}% lines, {tot.branch_pct:.1f}% branches")
    k = (m.result == "KILLED").sum(); s_ = (m.result == "SURVIVED").sum()
    b.bar(["killed", "survived"], [k, s_], color=[C["green"], C["red"]]); b.set_title(f"Mutation check: {k} of {k+s_} killed"); b.set_ylabel("hand-written mutants")
    for i, v in enumerate((k, s_)): b.text(i, v + 0.5, str(v), ha="center")
    if len(axs) == 3:
        t = pd.read_csv(R + "c8_test_counts.csv").sort_values("checks"); axs[2].barh(t.program, t.checks, color=C["dark"]); axs[2].set_xscale("log"); axs[2].set_xlabel("checks (plain run)"); axs[2].set_title(f"{int(t.checks.sum()):,} automated checks")
    src(f, "results/c8_coverage.csv, c8_mutation.csv" + (", c8_test_counts.csv" if len(axs) == 3 else "")); save(f, "V28", "verification")

# ---------------------------------------------------------------- G. optional
@fig("V29")
def v29():
    import json
    p = pd.read_csv(R + "c9_positions.csv"); meta = pd.read_csv(R + "c9_pair_meta.csv", header=None); A, B = int(meta.iloc[0, 1]), int(meta.iloc[0, 2])
    r = np.sqrt(p.x_km ** 2 + p.y_km ** 2 + p.z_km ** 2); p = p[(r > 6500) & (r < 8400)]
    pts = [[round(a, 1), round(b, 1), round(c, 1)] for a, b, c in zip(p.x_km, p.y_km, p.z_km)]
    hl = [[round(a, 1), round(b, 1), round(c, 1)] for a, b, c in zip(*[p[p.norad_id.isin([A, B])][k] for k in ("x_km", "y_km", "z_km")])]
    html = """<!doctype html><meta charset=utf-8><title>LEO objects at one conjunction</title><style>body{margin:0;background:#0b1220;color:#dde;font:13px sans-serif}canvas{display:block;margin:auto}#t{position:absolute;left:12px;top:8px}</style>
<div id=t>Drag to rotate, wheel to zoom. __N__ LEO objects at the closest crossing encounter of the 1 h run (SGP4, dlscms prop). Red = pair __A__ / __B__.</div><canvas id=c></canvas>
<script>const P=__P__,H=__H__;const cv=document.getElementById('c'),g=cv.getContext('2d');let W,Ht,ax=0.4,ay=0.6,z=1;function rs(){W=cv.width=innerWidth;Ht=cv.height=innerHeight}rs();onresize=rs;
let d=false,lx=0,ly=0;cv.onmousedown=e=>{d=true;lx=e.clientX;ly=e.clientY};onmouseup=()=>d=false;onmousemove=e=>{if(d){ay+=(e.clientX-lx)*.005;ax+=(e.clientY-ly)*.005;lx=e.clientX;ly=e.clientY;draw()}};onwheel=e=>{z*=e.deltaY<0?1.1:.9;draw();e.preventDefault()};
function pr(p){let x=p[0],y=p[1],zz=p[2],c=Math.cos(ay),s=Math.sin(ay);let x1=x*c-y*s,y1=x*s+y*c;c=Math.cos(ax);s=Math.sin(ax);let y2=y1*c-zz*s,z2=y1*s+zz*c;const k=Math.min(W,Ht)/18000*z;return[W/2+x1*k,Ht/2-z2*k,y2]}
function draw(){g.clearRect(0,0,W,Ht);const k=Math.min(W,Ht)/18000*z;g.fillStyle='#1d3b66';g.beginPath();g.arc(W/2,Ht/2,6378*k,0,7);g.fill();g.fillStyle='rgba(220,230,255,.55)';for(const p of P){const q=pr(p);if(q[2]>0||Math.hypot(q[0]-W/2,q[1]-Ht/2)>6378*k)g.fillRect(q[0],q[1],1.2,1.2)}
g.strokeStyle='#ff5a36';g.lineWidth=2;for(const p of H){const q=pr(p);g.beginPath();g.arc(q[0],q[1],7,0,7);g.stroke()}}draw()</script>"""
    html = html.replace("__P__", json.dumps(pts, separators=(",", ":"))).replace("__H__", json.dumps(hl)).replace("__N__", f"{len(p):,}").replace("__A__", str(A)).replace("__B__", str(B))
    open(OUT + "V29_orbit_viewer.html", "w").write(html); print("wrote", OUT + "V29_orbit_viewer.html")

@fig("V30")
def v30():
    if not os.path.exists(R + "c9_demo_log.txt"): print("V30 skipped: run scripts/demo.sh > results/c9_demo_log.txt first"); return
    lines = open(R + "c9_demo_log.txt").read().rstrip().split("\n")[:44]
    f = plt.figure(figsize=(11, 0.19 * len(lines) + 0.5)); f.patch.set_facecolor("#111418")
    f.text(0.01, 0.98, "\n".join(lines), family="monospace", fontsize=7.2, color="#d8dee9", va="top", ha="left")
    save(f, "V30", "demo_screenshot")

if __name__ == "__main__":
    want = sys.argv[1:] or sorted(FIGS)
    for k in want:
        try: FIGS[k]()
        except Exception as e: print("FAILED", k, type(e).__name__, e)
