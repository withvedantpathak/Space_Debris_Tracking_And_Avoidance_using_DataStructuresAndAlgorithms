"""C6 figures from results/c6_*.csv -> results/c6_summary.png"""
import os, pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R = lambda f: os.path.join(root, "results", f)
fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
# A: cluster size distribution (6 h)
ev = pd.read_csv(R("c6_events_6h.csv")); ev["a"] = ev[["norad_id_a", "norad_id_b"]].min(axis=1); ev["b"] = ev[["norad_id_a", "norad_id_b"]].max(axis=1)
pr = ev.groupby(["a", "b"]).size().reset_index()[["a", "b"]]
ids = pd.unique(pr.values.ravel()); ix = {v: i for i, v in enumerate(ids)}; par = list(range(len(ids)))
def find(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for a, b in zip(pr.a.map(ix), pr.b.map(ix)): par[find(a)] = find(b)
sizes = pd.Series([find(i) for i in range(len(ids))]).value_counts().values
vc = pd.Series(sizes).value_counts().sort_index()
ax[0].loglog(vc.index, vc.values, "o", color="#c2410c"); ax[0].set(xlabel="objects in cluster", ylabel="number of clusters", title=f"Clusters of conjunctions, 6 h window\n(largest: {sizes.max():,} objects)"); ax[0].grid(alpha=.3, which="both")
# B: cover sizes
s = pd.read_csv(R("c6_summary.csv")); s = s[s.variant == "reductions"].reset_index(drop=True)
lab = ["1 h", "3 h", "6 h"]; milp = [pd.read_csv(R(f"c6_{w}_verify.csv")).milp_optimum.iloc[0] for w in ["1h", "3h", "6h"]]
x = np.arange(3); wd = 0.2
for k, (name, vals, c) in enumerate([("matching lower bound", s.lower_bound, "#bbbbbb"), ("independent MILP optimum", milp, "#2563eb"), ("our exact search", s.exact_cover, "#1b7f5c"), ("greedy", s.greedy_cover, "#c2410c")]):
    ax[1].bar(x + (k - 1.5) * wd, vals, wd, label=name, color=c)
ax[1].set_xticks(x); ax[1].set_xticklabels(lab); ax[1].set(ylabel="objects to move (min vertex cover)", title="Mitigation: objects to manoeuvre"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, axis="y")
lo = min(s.lower_bound) * 0.85; ax[1].set_ylim(lo, max(s.greedy_cover) * 1.03)
# C: exact-search time per cluster
c = pd.concat([pd.read_csv(R(f"c6_{w}_components.csv")).assign(window=w) for w in ["1h", "3h", "6h"]])
c["t"] = c.exact_us.clip(lower=0.5)
for done, col, nm in [(1, "#1b7f5c", "search finished"), (0, "#c2410c", "cut off at node limit")]:
    q = c[c.exact_complete == done]; ax[2].scatter(q.nodes, q.t, s=14, color=col, alpha=.6, label=nm)
ax[2].set(xscale="log", yscale="log", xlabel="objects in cluster", ylabel="exact search time (microseconds)", title="Exact search cost per cluster (3+ objects)"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3, which="both")
fig.tight_layout(); fig.savefig(R("c6_summary.png"), dpi=130)
print(s[["events_file", "pairs", "nodes", "components", "largest_nodes", "greedy_cover", "exact_cover", "lower_bound"]].to_string(index=False)); print("MILP optimum:", milp)
