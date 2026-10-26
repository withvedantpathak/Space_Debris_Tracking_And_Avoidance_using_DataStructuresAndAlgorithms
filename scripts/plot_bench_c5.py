"""Plot results/bench_c5.csv -> results/bench_c5.png and results/bench_c5_summary.csv (speed-ups vs brute force)."""
import os, pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
d = pd.read_csv(os.path.join(root, "results", "bench_c5.csv"))
assert (d.identical_to_bruteforce == "yes").all(), "a fast method disagreed with brute force"
syn = d[d.dataset.str.contains("synth")]
col = {"brute": "#444444", "kd": "#c2410c", "grid": "#1b7f5c", "alt": "#2563eb"}
lab = {"brute": "brute force O(n^2)", "kd": "k-d tree per step", "grid": "hashed grid per step", "alt": "altitude filter + sweep"}
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for m in ["brute", "kd", "grid", "alt"]:
    s = syn[syn.method == m].sort_values("n")
    ax[0].loglog(s.n, s.cpu_seconds, "o-", color=col[m], label=lab[m]); ax[1].loglog(s.n, s.pair_checks, "o-", color=col[m], label=lab[m])
ax[0].set(xlabel="objects", ylabel="CPU seconds (1 h window, 10 s step, 10 km)", title="Time, synthetic catalogue"); ax[0].legend(); ax[0].grid(alpha=.3, which="both")
ax[1].set(xlabel="objects", ylabel="pairs examined", title="Pairs examined (log scale)"); ax[1].grid(alpha=.3, which="both")
r = d[d.dataset == "real_sgp4_1h"].set_index("method")
ms = [m for m in ["brute", "kd", "grid"] if m in r.index]
ax[2].bar([lab[m].split(" ")[0] + " " + lab[m].split(" ")[1] for m in ms], [r.loc[m, "cpu_seconds"] for m in ms], color=[col[m] for m in ms])
for i, m in enumerate(ms): ax[2].text(i, r.loc[m, "cpu_seconds"] + 1, f'{r.loc[m, "cpu_seconds"]:.0f} s', ha="center")
ax[2].set(ylabel="CPU seconds", title="Real catalogue, SGP4, 18,191 objects, 1 h, 10 km"); ax[2].grid(alpha=.3, axis="y")
fig.tight_layout(); fig.savefig(os.path.join(root, "results", "bench_c5.png"), dpi=130)
rows = []
for (ds, n), g in d.groupby(["dataset", "n"]):
    b = g[g.method == "brute"].cpu_seconds.iloc[0]
    for _, x in g.iterrows(): rows.append(dict(dataset=ds, n=n, method=x.method, cpu_seconds=x.cpu_seconds, speedup_vs_brute=round(b / x.cpu_seconds, 2), pair_checks=x.pair_checks, pair_check_reduction=round(g[g.method == "brute"].pair_checks.iloc[0] / x.pair_checks, 1)))
pd.DataFrame(rows).to_csv(os.path.join(root, "results", "bench_c5_summary.csv"), index=False)
print(pd.DataFrame(rows).to_string(index=False))
