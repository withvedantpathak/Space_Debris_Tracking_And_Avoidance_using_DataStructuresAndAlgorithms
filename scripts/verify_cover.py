"""Independent check of the C6 mitigation results with an integer-programming solver (scipy.optimize.milp, HiGHS).
For the conjunction graph in an events CSV, solve  min sum x_v  s.t.  x_u + x_v >= 1 for every clearable pair, x_v = 0 for immovable
objects, x binary, and compare with the C results (<prefix>_cover.csv): validity of both covers and the optimum size.
Usage: verify_cover.py <events.csv> <movable.csv|-> <prefix> -> prints a verdict and writes results/c6_verify.csv"""
import sys, os, time, numpy as np, pandas as pd
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
ev = pd.read_csv(sys.argv[1]); mov_path, prefix = sys.argv[2], sys.argv[3]
ev["a"] = ev[["norad_id_a", "norad_id_b"]].min(axis=1); ev["b"] = ev[["norad_id_a", "norad_id_b"]].max(axis=1)
pairs = ev.groupby(["a", "b"], as_index=False).miss_km.min()
ids = pd.unique(pairs[["a", "b"]].values.ravel()); idx = {v: i for i, v in enumerate(ids)}; n = len(ids)
if mov_path == "-": mov = np.ones(n, bool)
else:
    m_ids = set(pd.read_csv(mov_path, usecols=["NORAD_CAT_ID"]).NORAD_CAT_ID); mov = np.array([i in m_ids for i in ids])
u = pairs.a.map(idx).values; v = pairs.b.map(idx).values
clear = mov[u] | mov[v]; u, v = u[clear], v[clear]; m = len(u)
A = coo_matrix((np.ones(2 * m), (np.r_[np.arange(m), np.arange(m)], np.r_[u, v])), shape=(m, n)).tocsr()
t0 = time.time()
res = milp(c=np.ones(n), constraints=LinearConstraint(A, lb=1, ub=np.inf), integrality=np.ones(n), bounds=Bounds(0, np.where(mov, 1, 0)))
t = time.time() - t0
assert res.success, res.message
opt = int(round(res.fun))
cov = pd.read_csv(prefix + "_cover.csv"); out = []
for col in ["in_greedy_cover", "in_best_cover"]:
    sel = set(cov.norad_id[cov[col] == 1]); x = np.array([i in sel for i in ids])
    valid = bool((x[u] | x[v]).all() and not (x & ~mov).any())
    out.append((col, int(x.sum()), valid))
print(f"clearable pairs {m}, objects {n}, MILP optimum {opt} ({t:.2f} s, status: {res.message.strip()})")
for col, size, valid in out: print(f"  {col}: size {size}, valid cover: {valid}, gap to optimum: {size - opt}")
pd.DataFrame([dict(events=sys.argv[1], pairs=m, milp_optimum=opt, milp_seconds=round(t, 3), **{f"{c}_size": s for c, s, _ in out}, **{f"{c}_valid": v_ for c, _, v_ in out})]).to_csv(prefix + "_verify.csv", index=False)
comp = pd.read_csv(prefix + "_components.csv"); all_done = bool((comp.exact_complete == 1).all())
valid = out[1][2] and out[0][2] and out[0][1] >= opt and out[1][1] >= opt           # no cover can beat the true optimum
best_ok = valid and (out[1][1] == opt or not all_done)
print(f"  exact search finished on every cluster: {all_done}")
print("RESULT:", ("PASS (exact cover equals the MILP optimum)" if out[1][1] == opt else "PASS (search was cut off on a cluster; cover valid, optimum not proven, gap %d)" % (out[1][1] - opt)) if best_ok else "FAIL")
sys.exit(0 if best_ok else 1)
