"""Plot results/bench_c1.csv: runtime, speed-up, memory, recall (mean and min-max over seeds)."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src, out = sys.argv[1], sys.argv[2]
df = pd.read_csv(src)
df["total_ms"] = df["build_ms"] + df["query_ms"]
R = df["radius_km"].iloc[0]
col = {"brute": "#888888", "grid": "#1b7f5c", "kdtree": "#c2410c"}
g = df.groupby(["method", "n"]).agg(total=("total_ms", "mean"), tmin=("total_ms", "min"), tmax=("total_ms", "max"),
                                     build=("build_ms", "mean"), query=("query_ms", "mean"),
                                     mem=("mem_bytes", "mean"), recall=("recall", "min"), pairs=("pairs", "mean")).reset_index()
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
for m, d in g.groupby("method"):
    ax[0][0].errorbar(d["n"], d["total"], yerr=[d["total"] - d["tmin"], d["tmax"] - d["total"]], marker="o", label=m, color=col[m], capsize=3)
ax[0][0].set(xscale="log", yscale="log", xlabel="objects N", ylabel="build + all queries (ms)", title="Runtime vs N")
b = g[g.method == "brute"].set_index("n")["total"]
for m in ("grid", "kdtree"):
    d = g[g.method == m].set_index("n")
    ax[0][1].plot(d.index, b[d.index] / d["total"], marker="o", label=m, color=col[m])
ax[0][1].set(xscale="log", yscale="log", xlabel="objects N", ylabel="speed-up over brute force (x)", title="Speed-up")
for m in ("grid", "kdtree"):
    d = g[g.method == m]
    ax[1][0].plot(d["n"], d["mem"] / 1e6, marker="o", label=m, color=col[m])
ax[1][0].set(xscale="log", yscale="log", xlabel="objects N", ylabel="index memory (MB)", title="Memory")
for m in ("grid", "kdtree"):
    d = g[g.method == m]
    ax[1][1].plot(d["n"], d["recall"] * 100, marker="o", label=m, color=col[m])
ax[1][1].set(xscale="log", ylim=(95, 101), xlabel="objects N", ylabel="worst-case recall vs brute force (%)", title="Accuracy (min over seeds)")
for a in ax.flat:
    a.grid(alpha=.3); a.legend()
fig.suptitle(f"C1 micro-benchmark: pairs within {R:g} km, synthetic orbit-like shell, 3 seeds (mean, min-max)")
fig.tight_layout()
fig.savefig(out, dpi=130)
print(g.round(3).to_string(index=False))
